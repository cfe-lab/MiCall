from base64 import standard_b64encode
from io import BytesIO
from pathlib import Path
from turtle import Turtle

# noinspection PyPep8Naming
import drawsvg as draw
# noinspection PyPackageRequirements
from PIL import Image
from drawsvg import Drawing


def drawing_to_image(drawing: Drawing) -> Image:
    png = drawing.rasterize(context=draw.Context(invert_y=True))
    png_bytes = BytesIO(png.png_data)
    image = Image.open(png_bytes)
    return image


def encode_image(image: Image) -> str:
    writer = BytesIO()
    image.save(writer, format='PNG')
    encoded = standard_b64encode(writer.getvalue())
    return encoded.decode('UTF-8')


class SvgDiffer:
    def __init__(self):
        self.work_dir: Path = Path(__file__).parent / 'svg_diffs'
        self.work_dir.mkdir(exist_ok=True)
        self.mismatch_found = False
        for work_file in self.work_dir.iterdir():
            if work_file.name == '.gitignore':
                continue
            assert work_file.suffix in ('.svg', '.png')
            work_file.unlink()

    def diff_pixel(self, actual_pixel, expected_pixel):
        ar, ag, ab, aa = actual_pixel
        er, eg, eb, ea = expected_pixel
        if actual_pixel != expected_pixel:
            self.mismatch_found = True
            # Colour
            dr = 0xff
            dg = (ag + eg) // 5
            db = (ab + eb) // 5

            # Opacity
            da = 0xff
        else:
            # Colour
            dr, dg, db = ar, ag, ab

            # Opacity
            da = aa // 3
        return dr, dg, db, da

    def assert_equal(self,
                     svg_actual: Drawing,
                     svg_expected: Drawing,
                     name: str):
        png_actual = drawing_to_image(svg_actual)
        png_expected = drawing_to_image(svg_expected)
        w = max(png_actual.width, png_expected.width)
        h = max(png_actual.height, png_expected.height)

        png_actual_padded = Image.new(png_actual.mode, (w, h))
        png_expected_padded = Image.new(png_expected.mode, (w, h))
        png_actual_padded.paste(png_actual)
        png_expected_padded.paste(png_expected)
        png_diff = Image.new(png_actual.mode, (w, h))
        self.mismatch_found = False
        # noinspection PyTypeChecker
        png_diff.putdata([self.diff_pixel(actual_pixel, expected_pixel)
                          for actual_pixel, expected_pixel in zip(
                            png_actual_padded.get_flattened_data(),
                            png_expected_padded.get_flattened_data())])

        # Display image when in live turtle mode.
        display_image = getattr(Turtle, 'display_image', None)
        if display_image is not None:
            t = Turtle()
            try:
                # noinspection PyUnresolvedReferences
                w = t.screen.cv.cget('width')
                # noinspection PyUnresolvedReferences
                h = t.screen.cv.cget('height')
                ox, oy = w/2, h/2
                text_height = 20
                t.penup()
                t.goto(-ox, oy)
                t.right(90)
                t.forward(text_height)
                t.write(f'Actual')
                display_image(position=t.position(),
                              image=encode_image(png_actual))
                t.forward(png_actual.height)
                t.forward(text_height)
                t.write(f'Diff')
                display_image(position=t.position(),
                              image=encode_image(png_diff))
                t.forward(png_diff.height)
                t.forward(text_height)
                t.write('Expected')
                display_image(position=t.position(),
                              image=encode_image(png_expected))
                t.forward(png_expected.height)
            except Exception as ex:
                t.write(str(ex))

        if not self.mismatch_found:
            return
        text_actual = svg_actual.as_svg(context=draw.Context(invert_y=True))
        (self.work_dir / (name+'_actual.svg')).write_text(text_actual)
        text_expected = svg_expected.as_svg(context=draw.Context(invert_y=True))
        (self.work_dir / (name+'_expected.svg')).write_text(text_expected)
        with (self.work_dir / (name+'_diff.png')) as f:
            png_diff.save(f)
        assert text_actual == text_expected
