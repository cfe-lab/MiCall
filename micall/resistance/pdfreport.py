#!/usr/bin/env python3.6

# The module that generates a report in PDF format
import pytz
import datetime

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import cm, mm
from reportlab.lib.styles import ParagraphStyle
from reportlab import rl_config

import reportlab.platypus as plat

# we currently only support North American letter paper -- no A4
from reportlab.platypus import PageBreak

page_w, page_h = letter

# Make PDF's repeatable - don't generate id and force timestamp to 01 Jan 2000.
rl_config.invariant = 1

# Times are reported in this time zone
TIME_ZONE_NAME = "America/Vancouver"
time_zone = pytz.timezone(TIME_ZONE_NAME)

# text size in the table in points
TAB_FONT_SIZE = 8

# text size in the small print
SMALL_PRINT_FONT_SIZE = 7

# a style used for the 'relevant mutations' text
mut_txt_style = ParagraphStyle("sconormal",
                               fontSize=TAB_FONT_SIZE,
                               leading=TAB_FONT_SIZE,
                               fontName='Helvetica-Oblique')


def get_now_string():
    """Return the date and time in the configured time zone as a string"""
    utc_now = datetime.datetime.now(tz=pytz.utc)
    loc_now = utc_now.astimezone(time_zone)
    # return loc_now.strftime('%Y-%b-%d %H:%M:%S %Z')
    return loc_now.strftime('%Y-%b-%d (%Z)')


def bottom_para(txt):
    """Set the provided text into a form for the small print"""
    small_style = ParagraphStyle("small",
                                 fontSize=SMALL_PRINT_FONT_SIZE,
                                 leading=SMALL_PRINT_FONT_SIZE-1)
    return plat.Paragraph(txt, small_style)


def test_details_para(txt):
    """Set the provided text into a form for the test details"""
    small_style = ParagraphStyle("small",
                                 fontSize=TAB_FONT_SIZE,
                                 leading=TAB_FONT_SIZE-1)
    return plat.Paragraph(txt, small_style)


def headertab_style(row_offset, colnum, dospan):
    """Generate a style list for the first row of a table with colnum columns.
    dospan := turn the colnum columns into a single one with centred text

    This routine is responsible for the style in the table headings.
    """
    lst = [("TEXTCOLOR", (0, row_offset), (colnum-1, row_offset), colors.white),
           ("BACKGROUND", (0, row_offset), (colnum-1, row_offset), colors.green),
           ("ALIGN", (0, row_offset), (colnum-1, row_offset), "CENTRE"),
           ("FACE", (0, row_offset), (colnum-1, row_offset), "Helvetica-Bold")]
    if dospan:
        lst.extend([("SPAN", (0, row_offset), (colnum-1, row_offset))])
        # ("BOX", (0, row_offset), (colnum-1, row_offset), 1, colors.black)])
    else:
        lst.extend([("GRID", (0, row_offset), (colnum-1, row_offset), 0.5, colors.black)])
    return lst


def drug_class_tablst(row_offset, report_template, genotype, drug_class_code, level_coltab):
    cfg_dct = report_template.virus_config
    drug_lst = cfg_dct["known_drugs"][drug_class_code]
    table_header_str = cfg_dct['drug_class_tableheaders'][drug_class_code]
    report_page = report_template.genotype_pages[genotype]
    resistance_dct = report_page.resistance_calls
    mutation_str = report_page.mutations.get(drug_class_code, 'None')
    mutation_str = "Relevant {} Mutations: {}".format(drug_class_code,
                                                      mutation_str)
    # 1) row 0: header column: name of drug_class
    t_data = [["{} Drugs".format(table_header_str), ""]]
    t_style = headertab_style(row_offset, 2, dospan=True)
    # 2) row 1..num_drugs: list of drugs in this drug_class
    drow_min, drow_max = row_offset + 1,  row_offset + len(drug_lst)
    t_style.append(("GRID", (0, drow_min), (1, drow_max), 1.0, colors.white))
    t_style.extend([("ALIGNMENT", (0, drow_min), (0, drow_max), 'LEFT'),
                    ("LEFTPADDING", (0, drow_min), (0, drow_max), 0),
                    ("ALIGNMENT", (1, drow_min), (1, drow_max), 'CENTRE'),
                    ("FACE", (1, drow_min), (1, drow_max), "Helvetica-Bold")])
    for tabline, dd in enumerate(drug_lst):
        drug_id, drug_name = dd
        if drug_id in resistance_dct:
            level, level_name = resistance_dct[drug_id]
        else:
            level, level_name = 0, "Not indicated: genotype " + genotype
        t_data.append([drug_name, level_name])
        # determine colours for the level
        bg_col, fg_col = level_coltab[level]
        t_style.extend([('TEXTCOLOR', (1, tabline + drow_min), (1, tabline + drow_min), fg_col),
                        ('BACKGROUND', (1, tabline + drow_min), (1, tabline + drow_min), bg_col)])
        # if compact and tabline % 2 == 0:
        #    t_style.append(('BACKGROUND', (0, tabline + drow_min), (0, tabline + drow_min), colors.lightgrey))
    # 3) mutation string
    # we put this into a separate paragraph into a column that spans the two table columns
    mut_row = drow_max + 1
    t_style.extend([("SPAN", (0, mut_row), (1, mut_row)),
                    ("LEFTPADDING", (0, mut_row), (1, mut_row), 0)])
    # ("BOX", (0, mut_row), (1, mut_row), 0.5, colors.black)])
    t_data.append([plat.Paragraph(mutation_str, mut_txt_style), ""])
    assert sum([len(row) == 2 for row in t_data]) == len(t_data), "wonky drug table"
    return t_data, t_style


def top_table(sample_name, table_width, genotype):
    """Generate a (mostly empty) top table of three main columns.
    table_width: the overall width of the table.
    """
    samp_name = sample_name or "None"
    mid_colwidth = table_width/2.8
    oth_colwidth = (table_width - mid_colwidth)/2.0
    nowstr = get_now_string()
    test_dl = [["Patient/Sample Details", "Test Details", "Physician Details"],
               ["", test_details_para("Sample ID: {}".format(samp_name)), ""],
               ["", test_details_para("Report Date: {}".format(nowstr)), ""],
               ["",
                (genotype or "") and
                test_details_para("Genotype: " + genotype),
                ""],
               ["", "", ""]
               ]
    rn_min, rn_max = 1, len(test_dl) - 1
    lc, mc, rc = 0, 1, 2
    st_lst = headertab_style(0, 3, dospan=False)
    st_lst.extend([("BOX", (lc, rn_min), (lc, rn_max), 0.5, colors.black),
                   ("BOX", (rc, rn_min), (rc, rn_max), 0.5, colors.black),
                   ("GRID", (mc, rn_max+1), (mc, rn_max), 0.5, colors.black),
                   ("FONTSIZE", (lc, rn_min), (rc, rn_max), 8)])
    return plat.Table(test_dl, style=st_lst,
                      colWidths=[oth_colwidth, mid_colwidth, oth_colwidth],
                      hAlign="CENTRE")


def write_report_one_column(report_templates, fname, sample_name=None):
    """Generate a PDF report to a given output file name
    """
    doc = plat.SimpleDocTemplate(
        fname,
        pagesize=letter,
        topMargin=1 * cm,
        title="basespace drug resistance report",
        author="BCCfE")
    # get the actual text width, (not the page width):
    # noinspection PyUnresolvedReferences
    txt_w = page_w - doc.leftMargin - doc.rightMargin
    table_width = txt_w - 1 * cm
    doc_els = []
    ti_style = ParagraphStyle("scotitle", alignment=TA_CENTER, fontSize=20)
    re_style = ParagraphStyle("scored", fontSize=15, textColor=colors.red,
                              spaceBefore=5 * mm, spaceAfter=5 * mm)
    failure_template = None
    for report_template in report_templates:
        if 'failure_message' in report_template.virus_config:
            failure_template = report_template
            continue
        # from the resistance, we determine which drug_classes to write a table for:
        # we only write a table if we are given resistance data for it.
        reported_genotypes = report_template.get_reported_genotypes()
        if not reported_genotypes:
            continue
        for genotype in reported_genotypes:
            if doc_els:
                doc_els.append(PageBreak())
            got_dc_set = report_template.get_reported_drug_classes(genotype)
            cfg_dct = report_template.virus_config
            col_tab = cfg_dct["resistance_level_colours"]
            level_coltab = dict([(k, (colors.HexColor(v[1]), colors.HexColor(v[2])))
                                 for k, v in col_tab.items()])
            append_header(doc_els,
                          cfg_dct,
                          table_width,
                          re_style,
                          ti_style,
                          sample_name,
                          genotype)
            # now drug classes tables, two per line
            known_dc_lst = cfg_dct["known_dclass_list"]
            tot_tab, tot_style = [], []
            for dc in [dc for dc in known_dc_lst if dc in got_dc_set]:
                tl, t_style = drug_class_tablst(len(tot_tab),
                                                report_template,
                                                genotype,
                                                dc,
                                                level_coltab)
                tot_tab.extend(tl)
                tot_style.extend(t_style)
            # adjust the overall table style
            num_rows = len(tot_tab)
            tot_style.extend([("VALIGN", (0, 0), (1, num_rows-1), "TOP"),
                              ("FONTSIZE", (0, 0), (1, num_rows-1), TAB_FONT_SIZE),
                              ("LEADING", (0, 0), (1, num_rows-1), TAB_FONT_SIZE)])
            left_col_w = table_width * 0.5
            right_col_w = table_width - left_col_w
            doc_els.append(plat.Table(tot_tab,
                                      vAlign="TOP",
                                      hAlign="CENTRE", style=tot_style,
                                      colWidths=[left_col_w, right_col_w]))
            # this is for layout debugging
            # big_table = [["l0", "r0"], ["l1", "r1"], ["l2", "r2"]]
            # debug_lst = [("GRID", (lc, 0), (rc, d_rowmax), 1, colors.red)]
            # btstyle.extend(debug_lst)
            append_footer(doc_els, cfg_dct)
    assert failure_template is not None
    if not doc_els:
        cfg_dct = failure_template.virus_config
        append_header(doc_els,
                      cfg_dct,
                      table_width,
                      re_style,
                      ti_style,
                      sample_name)
        doc_els.append(plat.Paragraph(
            "Sequence does not meet quality-control standards",
            re_style))
        append_footer(doc_els, cfg_dct)
    doc.build(doc_els)


def append_footer(doc_els, cfg_dct):
    # bottom paragraphs
    for field_name in ('disclaimer_text', 'generated_by_text'):
        field_text = cfg_dct.get(field_name)
        if field_text:
            doc_els.append(bottom_para(field_text))


def append_header(doc_els,
                  cfg_dct,
                  table_width,
                  re_style,
                  ti_style,
                  sample_name,
                  genotype=None):
    doc_els.append(plat.Paragraph(cfg_dct["report_title"], ti_style))
    doc_els.append(plat.Paragraph("For research use only", re_style))
    # -- top table
    doc_els.append(top_table(sample_name, table_width, genotype))


if __name__ == '__main__':
    print("The local time is '{}'".format(get_now_string()))
    # gen_testpage("testpage.pdf")
