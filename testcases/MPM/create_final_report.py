from pypdf import PdfReader, PdfWriter, PageObject, PdfWriter
import img2pdf
import os

#-------------------------------------------------------------------------------

def stitch_pdfs_side_by_side(pdf1_path,
                             pdf2_path,
                             output_path):

    reader1 = PdfReader(pdf1_path)
    reader2 = PdfReader(pdf2_path)
    writer = PdfWriter()

    # Grab the first page from each PDF
    page1 = reader1.pages[0]
    page2 = reader2.pages[0]

    # Extract dimensions
    width1 = page1.mediabox.width
    height1 = page1.mediabox.height
    width2 = page2.mediabox.width
    height2 = page2.mediabox.height

    # Calculate dimensions for the wide canvas
    total_width = width1 + width2
    max_height = max(height1, height2)

    # Create the modern blank layout page
    new_page = PageObject.create_blank_page(width=total_width, height=max_height)

    # Merge pages onto the canvas
    new_page.merge_page(page1)
    new_page.merge_translated_page(page2, tx=width1, ty=0)

    # Write out the final canvas
    writer.add_page(new_page)
    with open(output_path, "wb") as f:
        writer.write(f)

#-------------------------------------------------------------------------------

def create_final_report():

    """
    This function combines testcase plots from this checkout of the testcases with
    a baseline checkout of the testcases
    """

    testcases = {"polynya1":{"testcaseLocation":"polynya",
                             "plotFilename":"cells_scatter.pdf"},
                 "polynya2":{"testcaseLocation":"polynya",
                             "plotFilename":"particles_scatter.pdf"},
                 "interpolation":{"testcaseLocation":"spherical_operators/interpolation",
                                  "plotFilename":"interpolation_scaling.pdf"},
                 "reconstruction1":{"testcaseLocation":"spherical_operators/reconstruction",
                                    "plotFilename":"reconstruction_scaling.pdf"},
                 "reconstruction2":{"testcaseLocation":"spherical_operators/reconstruction",
                                    "plotFilename":"reconstruction_map.pdf"},
                 "interpolation_cell":{"testcaseLocation":"spherical_operators/interpolation_cell",
                                       "plotFilename":"interpolation_cell_scaling.pdf"}}

    writer = PdfWriter()

    for testcase in testcases.keys():

        filename1 = "./" + testcases[testcase]["testcaseLocation"] + "/" + testcases[testcase]["plotFilename"]
        filename2 =  "./baseline_plots/" + testcases[testcase]["plotFilename"]

        stitch_pdfs_side_by_side(filename1, filename2, "tmp.pdf")

        writer.append("tmp.pdf")

    with open("final_report.pdf", "wb") as f:
        writer.write(f)

    writer.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_final_report()
