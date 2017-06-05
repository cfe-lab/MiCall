#!/usr/bin/env python3.4

# The module that generates a report in PDF format
import datetime

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import cm
from reportlab.lib.styles import ParagraphStyle

import reportlab.platypus as plat

# stylesheet = getSampleStyleSheet()
# print("BLA", stylesheet.keys())
# normal_style = stylesheet["Normal"]
# small_style = stylesheet["Small"]

title_str = "HIV Protease/RT Resistance Genotype Report"

page_w, page_h = letter

top_margin = 2 * cm
logo_height = 2 * cm
title_pos = page_h - (top_margin + logo_height)


def wr_first_page(canv, doc):
    canv.saveState()
    canv.setFont('Times-Bold', 20)
    canv.drawCentredString(page_w / 2, title_pos, title_str)


def bottom_para(txt):
    small_style = ParagraphStyle("small", fontSize=8)
    return plat.Paragraph(txt, small_style)


def drug_class_table(cfg_dct, dc_name, level_coltab):
    """Generate a resistance report for a given drug class
    lr_align: str 'LEFT' or 'RIGHT': alignment of table on page.
    """
    drug_lst = cfg_dct["known_drugs"][dc_name]
    resistance_dct = cfg_dct["res_results"]
    mutation_str = cfg_dct["mutation_strings"][dc_name]
    # header column: name of drug_class
    t_data = [["{} Drugs".format(dc_name)]]
    num_drugs = len(drug_lst)
    t_style = [("BOX", (0, 0), (1, 0), 1, colors.black),
               ("GRID", (0, 1), (1, num_drugs), 1, colors.black),
               ("BOX", (0, num_drugs + 1), (1, num_drugs + 2), 1,
                colors.black),
               ("SPAN", (0, num_drugs + 1), (1, num_drugs + 1))]
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
    t_data.append(["Relevant {} Mutations:".format(dc_name), ""])
    t_data.append([mutation_str])
    t = plat.Table(t_data, vAlign="TOP")  # , hAlign=lr_align)
    t.setStyle(t_style)
    return t


def top_table(sample_name=None):
    """Generate a mostly empty top table."""
    samp_name = sample_name or "None"
    test_dl = [["Test Details"], ["Sample ID:", samp_name],
               ["Report Date (UTC):", str(datetime.datetime.utcnow())]]
    T_lines = len(test_dl)
    test_det_tab = plat.Table(test_dl)
    test_det_tab.setStyle([("BOX", (0, 0), (1, 0), 1, colors.black),
                           ("GRID", (0, 1), (1, T_lines), 1, colors.black)])
    bt = [["", test_det_tab, ""]]
    big_tab = plat.Table(bt)
    big_tab.setStyle([("GRID", (0, 0), (2, 0), 1, colors.black)])
    return big_tab


def write_report(cfg_dct, res_lst, mut_lst, fname, sample_name=None):
    """Generate a PDF report to a given output filename
    """
    col_tab = cfg_dct["resistance_level_colours"]
    level_coltab = dict([(k, (colors.HexColor(v[1]), colors.HexColor(v[2])))
                         for k, v in col_tab.items()])
    doc = plat.SimpleDocTemplate(
        fname,
        pagesize=letter,
        title="basespace drug resistance report",
        author="BCCfE")
    doc_els = []
    doc_els.append(plat.Spacer(1, 1.5 * cm))

    # a table that contains all other tables
    big_table = [["For research use only"]]
    # top table (mostly empty)
    big_table.append([top_table(sample_name)])
    # now drug classes tables, two per line
    tl = [
        drug_class_table(cfg_dct, dc, level_coltab)
        for dc in ["NRTI", "NNRTI", "PI", "INSTI"]
    ]
    big_table.append([tl[0], tl[1]])
    big_table.append([tl[2], tl[3]])

    bt = plat.Table(big_table)
    bt.setStyle([
        ('VALIGN', (0, 0), (1, 3), 'TOP'),
        ('SPAN', (0, 1), (1, 1)),
        ('SPAN', (0, 0), (1, 0)),
    ])
    doc_els.append(bt)
    doc_els.append(bottom_para(cfg_dct["disclaimer_text"]))
    small_print = cfg_dct["generated_by_text"] + "\n" + cfg_dct["version_text"]
    doc_els.append(bottom_para(small_print))
    doc.build(doc_els, onFirstPage=wr_first_page)


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
