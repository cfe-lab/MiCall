#!/usr/bin/env python3.4

# The module that generates a report in PDF format
import datetime

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import cm, mm
from reportlab.lib.styles import ParagraphStyle

import reportlab.platypus as plat

page_w, page_h = letter

top_margin = 2 * cm
logo_height = 2 * cm
title_pos = page_h - (top_margin + logo_height)


def bottom_para(txt):
    small_style = ParagraphStyle("small", fontSize=5, leading=6)
    return plat.Paragraph(txt, small_style)


def headertab_style(colnum, dospan):
    """Generate a style list for the first row of a table with colnum columns.
    dospan := turn the colnum columns into a single one with centred text"""
    lst = [
        ("TEXTCOLOR", (0, 0), (colnum-1, 0), colors.white),
        ("BACKGROUND", (0, 0), (colnum-1, 0), colors.green),
        ("GRID", (0, 0), (colnum-1, 0), 1, colors.black),
        ("ALIGN", (0, 0), (colnum-1, 0), "CENTRE")]
    if dospan:
        lst.extend([("SPAN", (0, 0), (colnum-1, 0)),
                    ("BOX", (0, 0), (colnum-1, 0), 1, colors.black)])
    return lst


def drug_class_table(cfg_dct, dc_name, level_coltab):
    """Generate a resistance report for a given drug class
    lr_align: str 'LEFT' or 'RIGHT': alignment of table on page.
    """
    drug_lst = cfg_dct["known_drugs"][dc_name]
    table_header_str = cfg_dct['drug_class_tableheaders'][dc_name]
    resistance_dct = cfg_dct["res_results"]
    mutation_str = cfg_dct["mutation_strings"][dc_name]
    # header column: name of drug_class
    num_drugs = len(drug_lst)
    drow_min, drow_max = 1,  num_drugs
    t_data = [["{} Drugs".format(table_header_str), ""]]
    t_style = [("GRID", (0, drow_min), (1, drow_max), 1, colors.black)]
    # header box style
    t_style.extend(headertab_style(2, dospan=True))
    # list of drugs in this drug_class
    for tabline, dd in enumerate(drug_lst):
        drug_id, drug_name = dd
        if drug_id in resistance_dct:
            level, level_name = resistance_dct[drug_id]
        else:
            level, level_name = 1, "NOT REPORTED"
        # determine colours for the level
        bg_col, fg_col = level_coltab[level]
        t_style.append(('TEXTCOLOR', (1, tabline + 1), (1, tabline + 1),
                        fg_col))
        t_style.append(('BACKGROUND', (1, tabline + 1), (1, tabline + 1),
                        bg_col))
        t_data.append([drug_name.capitalize(), level_name])
    # mutation string
    # a separate style list for the 'relevant mutations' box
    mut_row = drow_max + 1
    rel_mut_style = [("SPAN", (0, mut_row), (1, mut_row)),
                     ("BOX", (0, mut_row), (1, mut_row), 1, colors.black)]
    t_style.extend(rel_mut_style)
    mut_para_style = ParagraphStyle("sconormal")
    t_data.append([plat.Paragraph(mutation_str, mut_para_style), ""])
    assert sum([len(row) == 2 for row in t_data]) == len(t_data), "wonky drug table"
    return plat.Table(t_data, vAlign="TOP", style=t_style)  # , hAlign=lr_align)


def top_table(sample_name, colwidth):
    """Generate a mostly empty top table."""
    samp_name = sample_name or "None"
    # get the time, ignoring microseconds
    nowstr = str(datetime.datetime.utcnow().replace(microsecond=0))
    test_dl = [["Patient/Sample Details", "Test Details", "Physician Details"],
               ["", "Sample ID: {}".format(samp_name), ""],
               ["", "Report Date (UTC): {}".format(nowstr), ""],
               ["", "", ""]
               ]
    T_lines = len(test_dl)
    rn_min, rn_max = 1, T_lines - 1
    lc, mc, rc = 0, 1, 2
    st_lst = [("BOX", (lc, rn_min), (lc, rn_max), 1, colors.black),
              ("BOX", (rc, rn_min), (rc, rn_max), 1, colors.black),
              ("GRID", (mc, rn_min), (mc, rn_max), 1, colors.black),
              ("FONTSIZE", (0, rn_min), (2, rn_max), 8)]
    # colwidth -= 1.0 * cm
    st_lst.extend(headertab_style(3, dospan=False))
    tt = plat.Table(test_dl, style=st_lst, colWidths=[colwidth, None, colwidth])
    tt.hAlign = "CENTRE"
    return tt


def write_report(cfg_dct, res_lst, mut_lst, fname, sample_name=None):
    """Generate a PDF report to a given output file name
    """
    col_tab = cfg_dct["resistance_level_colours"]
    level_coltab = dict([(k, (colors.HexColor(v[1]), colors.HexColor(v[2])))
                         for k, v in col_tab.items()])
    doc = plat.SimpleDocTemplate(
        fname,
        pagesize=letter,
        title="basespace drug resistance report",
        author="BCCfE")
    # get the actual text width, (not the page width):
    txt_w = page_w - doc.leftMargin - doc.rightMargin
    w_half, w_top = txt_w * 0.5, txt_w / 3.0
    doc_els = [plat.Spacer(1, 1.5 * cm)]
    ti_style = ParagraphStyle("scotitle", alignment=TA_CENTER, fontSize=20)
    doc_els.append(plat.Paragraph(cfg_dct["report_title"], ti_style))
    re_style = ParagraphStyle("scored", fontSize=15, textColor=colors.red,
                              spaceBefore=5 * mm, spaceAfter=5 * mm)
    doc_els.append(plat.Paragraph("For research use only", re_style))
    # -- top table
    doc_els.append(top_table(sample_name, w_top))
    lc, rc = 0, 1
    big_table, btstyle = [], []
    # now drug classes tables, two per line
    known_dc_lst = cfg_dct["known_dclass_list"]
    # from the resistance, we determine which drug_classes to write a table for:
    # we only write a table if we are given resistance data for it.
    got_dc_set = set([dc["drug_class"] for dc in res_lst])
    tl = [drug_class_table(cfg_dct, dc, level_coltab) for dc in known_dc_lst if dc in got_dc_set]
    d_rowmin = 0
    while len(tl) > 0:
        row_lst = [tl.pop(0)]
        if len(tl) > 0:
            row_lst.append(tl.pop(0))
        else:
            row_lst.append("")
        big_table.append(row_lst)
    d_rowmax = len(big_table) - 1
    btstyle.extend([
        ("ALIGN", (lc, d_rowmin), (lc, d_rowmax), "RIGHT"),
        ("ALIGN", (rc, d_rowmin), (rc, d_rowmax), "LEFT"),
        ('VALIGN', (lc, d_rowmin), (rc, d_rowmax), 'TOP')])
    # this is for layout debugging
    # big_table = [["l0", "r0"], ["l1", "r1"], ["l2", "r2"]]
    # debug_lst = [("GRID", (lc, 0), (rc, d_rowmax), 1, colors.red)]
    # btstyle.extend(debug_lst)
    assert sum(len(row) == 2 for row in big_table) == len(big_table), "big_table lines are wonky"
    doc_els.append(plat.Table(big_table,
                              style=btstyle,
                              colWidths=[w_half, w_half],
                              hAlign="CENTRE"))
    # bottom paragraphs
    doc_els.append(bottom_para(cfg_dct["disclaimer_text"]))
    doc_els.append(bottom_para(cfg_dct["generated_by_text"]))
    doc.build(doc_els)


def gen_testpage(fname):
    write_report({}, [], [], fname)


def simple_gen_testpage(fname):
    """Generate a simple test page"""
    # NOTE: this example taken from
    # https://www.blog.pythonlibrary.org/2010/09/21/reportlab-tables-creating-tables-in-pdfs-with-python/
    doc = plat.SimpleDocTemplate(fname, pagesize=letter)
    elements = []
    data = [['00', '01', '02', '03', '04'], ['10', '11', '12', '13', '14'],
            ['20', '21', '22', '23', '24'], ['30', '31', '32', '33', '34']]
    t = plat.Table(data)
    t.setStyle(
        plat.TableStyle([('BACKGROUND', (1, 1), (-2, -2), colors.green), (
            'TEXTCOLOR', (0, 0), (1, -1), colors.red)]))
    elements.append(t)
    # write the document to disk
    doc.build(elements)


if __name__ == '__main__':
    gen_testpage("testpage.pdf")
