#!/usr/bin/env python3.4

# The module that generates a report in PDF format

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


def drug_class_table(cfg_dct, dc_name):
    """ Generate a resistance report for a given drug class"""
    t_data = [[dc_name]]
    drug_lst = cfg_dct["known_drugs"][dc_name]
    resistance_dct = cfg_dct["res_results"]
    for drug_id, drug_name in drug_lst:
        if drug_id in resistance_dct:
            level, level_name = resistance_dct[drug_id]
        else:
            level, level_name = 1, "NOT REPORTED"
        t_data.append([drug_name.capitalize(), level_name])
    t = plat.Table(t_data)
    return t


def write_report(cfg_dct, res_lst, mut_lst, fname):
    """Generate a PDF report to a given output filename
    """
    doc = plat.SimpleDocTemplate(
        fname,
        pagesize=letter,
        title="basespace drug resistance report",
        author="BCCfE")
    doc_els = []
    doc_els.append(plat.Spacer(1, 2 * cm))

    doc_els.append(drug_class_table(cfg_dct, 'NRTI'))
    doc_els.append(bottom_para(cfg_dct["disclaimer_text"]))
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
