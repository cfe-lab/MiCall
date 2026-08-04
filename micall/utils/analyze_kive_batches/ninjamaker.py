
from dataclasses import dataclass
from typing import Sequence, Union, NoReturn, Optional, Tuple
from pathlib import Path
import shlex


#
# Simple escaping for Ninja files.
# See <https://ninja-build.org/manual.html#ref_lexer> for explanation.
#
def escape(s: str) -> str:
    return s.replace('$', '$$').replace(' ', '$ ').replace('\n', '$\n').replace(':', '$:')


Value = Union[str, Path, 'Deref']
CommandArg = Union[str, Path, 'Deref']
CommandArgs = Union[CommandArg, Sequence[CommandArg]]


def compile_value(v: Value) -> str:
    if isinstance(v, Deref):
        return v.compile()
    elif isinstance(v, Path):
        return escape(str(v))
    elif isinstance(v, str):
        return escape(v)
    else:
        # Should never happen
        _: NoReturn = v
        raise RuntimeError(f"Unknown value type: {v!r}")


def compile_command_arg(a: CommandArg) -> str:
    """
    Use shlex.quote for shell safety, but preserve Deref.
    """
    if isinstance(a, Deref):
        return a.compile()
    elif isinstance(a, Path):
        return shlex.quote(escape(str(a)))
    elif isinstance(a, str):
        return shlex.quote(escape(a))
    else:
        _: NoReturn = a
        raise RuntimeError(f"Unexpected CommandArg type: {type(a)}.")


@dataclass(frozen=True)
class Deref:
    """
    A reference to a Ninja variable, e.g. ${out_dir}
    """
    name: str

    def compile(self) -> str:
        return f"${{{escape(self.name)}}}"


#
# == VARIABLE DEFINITIONS ==
#
@dataclass(frozen=True)
class VariableDefinition:
    name: str
    value: Value

    def compile(self) -> str:
        return f"{escape(self.name)} = {compile_value(self.value)}"


#
# -- COMMAND FOR A RULE --
#
@dataclass(frozen=True)
class Command:
    """
    head  := the program or ${var} to run
    arguments := a sequence of string|Path|Deref
    This will be rendered as:
      command = <head> <arg1> <arg2> …
    """

    head: Union[str, Deref]
    arguments: Sequence[CommandArg]

    def compile(self) -> str:
        head_s = compile_value(self.head)
        args_s = " ".join(compile_command_arg(a) for a in self.arguments)
        return head_s + ((" " + args_s) if args_s else "")

    @staticmethod
    def make(head: Union[str, Deref], *arguments: CommandArgs) -> 'Command':
        flattened: list[CommandArg] = []
        for arg in arguments:
            if isinstance(arg, (Deref, Path, str)):
                flattened.append(arg)
            else:
                flattened.extend(arg)

        return Command(head, tuple(flattened))


#
# -- DESCRIPTION FOR A RULE --
#
@dataclass(frozen=True)
class Description:
    format: str
    arguments: Sequence[Value]

    def compile(self) -> str:
        escaped_arguments = map(compile_value, self.arguments)
        return self.format.format(*escaped_arguments)

    @staticmethod
    def make(format: str, *arguments: Value) -> 'Description':
        return Description(format=format, arguments=arguments)


#
# == RULE DEFINITIONS ==
#
@dataclass(frozen=True)
class Rule:
    """
    rule <name>
      command = <cmd>
      description = <desc>   # optional
    """

    name: str
    command: Command
    description: Optional[Description] = None

    def compile(self) -> str:
        lines = [f"rule {self.name}", f"  command = {self.command.compile()}"]
        if self.description:
            lines.append(f"  description = {self.description.compile()}")
        return "\n".join(lines)


#
# == BUILD STATEMENTS ==
#
@dataclass(frozen=True)
class Build:
    """
    build <outs>: <rulename> <ins> | <implicit> || <order-only>
      var1 = value1
      var2 = value2
      ...
      varn = valuen
    """

    outputs: Sequence[Value]
    rule: str
    inputs: Sequence[Value]
    bindings: Sequence[Tuple[str, Value]] = ()
    implicit: Sequence[Value]             = ()
    order_only: Sequence[Value]           = ()

    def compile(self) -> str:
        out_s = " ".join(compile_value(v) for v in self.outputs)
        in_s  = " ".join(compile_value(v) for v in self.inputs)
        parts = [f"build {out_s}: {self.rule}"]
        if in_s:
            parts.append(in_s)
        if self.implicit:
            parts.append("|")
            parts.extend(compile_value(v) for v in self.implicit)
        if self.order_only:
            parts.append("||")
            parts.extend(compile_value(v) for v in self.order_only)

        parts_s = " ".join(parts)
        whole = parts_s

        if self.bindings:
            bindings_s = "\n  ".join(VariableDefinition(name=name, value=value).compile()
                                     for name, value in self.bindings)
            whole += "\n  " + bindings_s

        return whole


#
# == DEFAULT TARGETS ==
#
@dataclass(frozen=True)
class Default:
    """
    default <outs...>
    """
    outputs: Sequence[Value]

    def compile(self) -> str:
        outs = " ".join(compile_value(v) for v in self.outputs)
        return f"default {outs}"


#
# == THE WHOLE NINJA FILE ==
#
Statement = Union[VariableDefinition, Rule, Build]


@dataclass(frozen=True)
class Recipe:
    statements: Sequence[Statement]
    default: Sequence[Value] = ()

    def compile(self) -> str:
        # join with a blank line between top‐level statements
        chunks = list(stmt.compile() for stmt in self.statements)
        if self.default:
            chunks.append(Default(self.default).compile())

        return "\n\n".join(chunks) + "\n"
