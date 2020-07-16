""" Spring beads align along a wire, elastically attracted to their targets.

Used when aligning sections of a contig to a reference sequence.
"""
import typing


class Bead:
    def __init__(self,
                 start: int,
                 end: int,
                 target_start: int = None,
                 target_end: int = None,
                 alignment=None,
                 skipped=0):
        self.start = self.display_start = start
        self.end = self.display_end = end
        self.target_start = target_start
        self.target_end = target_end
        self.alignment = alignment
        self.skipped = skipped
        self.force = 0.0
        self.left: typing.Optional['Bead'] = None
        self.right: typing.Optional['Bead'] = None

    def __repr__(self):
        return f'Bead({self.start}, {self.end})'

    def set_target_force(self, k: float):
        if self.target_start is not None:
            bead_mid = (self.display_start + self.display_end) / 2
            target_mid = (self.target_start + self.target_end) / 2
            self.force = k*(target_mid - bead_mid)
        else:
            self.force = 0.0
            if self.right is not None:
                self.force += k*(self.right.display_start - self.display_end)
            if self.left is not None:
                self.force += k*(self.left.display_end - self.display_start)

    def apply_force(self) -> float:
        move = self.force
        left = self.left
        if left is not None and left.display_end == self.display_start:
            move += left.force
        right = self.right
        if right is not None and right.display_start == self.display_end:
            move += right.force
        if left is not None:
            min_move = left.display_end - self.display_start
            move = max(move, min_move)
        if right is not None:
            max_move = right.display_start - self.display_end
            move = min(move, max_move)
        self.display_start += move
        self.display_end += move
        return move


class Wire(list):
    def __repr__(self):
        beads = super().__repr__()
        return f'Wire({beads})'

    def align(self):
        prev_bead = None
        for bead in self:
            if prev_bead is not None:
                prev_bead.right = bead
                bead.left = prev_bead
            prev_bead = bead
        k = 0.5
        min_movement = k / 2
        while True:
            for bead in self:
                bead.set_target_force(k)
            movement = sum(abs(bead.apply_force()) for bead in self)
            if movement < min_movement:
                break
        for bead in self:
            bead.display_start = round(bead.display_start)
            bead.display_end = round(bead.display_end)

    @property
    def displays(self):
        return [(bead.display_start, bead.display_end)
                for bead in self]
