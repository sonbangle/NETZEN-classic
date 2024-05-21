import argparse,sys
#### import the simple module from the paraview
from paraview.simple import *
import os
import os.path


def export_network_image(infile='/data/SonData/BlenderExercises/vtk_test/data/layout/paraview_1.4.4/1.4.4_edgeTrue_subnet_and_gene.vtu',
                 outfile='/data/SonData/BlenderExercises/vtk_test/data/layout/paraview_1.4.4/1.4.4_edgeTrue_subnet_and_gene.png',
                 ImageResolution=[13990, 8380],
                 coloring_field='enrichment_subnet_score',
                 SetScalarBarVisibility=True,
                 min_val=0.0,
                 max_val=1.0,
                 opacity=0.7,
                 delete_input_file=True): # delete vtu input file after image generation

    #opacity =1 # always transparent to see the font
    #### disable automatic camera reset on 'Show'

    # Check if the outfile exists or not. if exists, then do nothing.
    if os.path.exists(outfile):
        print("output file" + outfile + " exists. Exit ")
        return None
    else:
        print("output file" + outfile + " does not exist. Generating image")
    paraview.simple._DisableFirstRenderCameraReset()
    # create a new 'XML Unstructured Grid Reader'
    data = XMLUnstructuredGridReader(FileName=[infile])
    data.PointArrayStatus = [coloring_field]
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    dataDisplay = Show(data, renderView1)
    # trace defaults for the display properties.
    dataDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera()

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(dataDisplay, ('POINTS', coloring_field))

    # rescale color and/or opacity maps used to include current data range
    dataDisplay.RescaleTransferFunctionToDataRange(True, False)



    # get color transfer function/color map for 'enrichment_subnet_score'
    LUT = GetColorTransferFunction(coloring_field)

    # Rescale transfer function
    LUT.RescaleTransferFunction(min_val, max_val)

    # get opacity transfer function/opacity map for 'enrichment_subnet_score'
    PWF = GetOpacityTransferFunction(coloring_field)

    # Rescale transfer function
    PWF.RescaleTransferFunction(min_val, max_val)

    # Properties modified on renderView1
    renderView1.Background = [0.0, 0.0, 0.0]

    dataDisplay.Opacity = opacity

    # get color legend/bar for enrichment_subnet_scoreLUT in view renderView1
    LUTColorBar = GetScalarBar(LUT, renderView1)

    # Properties modified on enrichment_subnet_scoreLUTColorBar
    LUTColorBar.TitleFontSize = 6
    LUTColorBar.LabelFontSize = 6
    LUTColorBar.ScalarBarThickness = 6

    # Properties modified on enrichment_subnet_scoreLUTColorBar
    LUTColorBar.ScalarBarLength = 0.2

    # show color bar/color legend
    dataDisplay.SetScalarBarVisibility(renderView1, SetScalarBarVisibility)

    # save screenshot
    SaveScreenshot(outfile,
        renderView1, ImageResolution=ImageResolution,
        ImageQuality=80)
    if delete_input_file:
        os.system("rm -f " + infile)
        os.system("rm -f " + infile + "_filter_condition.pickle")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--infile",help="input  vtu file", default='/data/SonData/BlenderExercises/vtk_test/data/layout/view_1.4.4/1.4.4_edgeFalse_gene_only1.vtu')
    parser.add_argument("--outfile",help="output png file", default='/data/SonData/BlenderExercises/vtk_test/data/layout/view_1.4.4/1.4.4_edgeFalse_gene_only1.png')
    parser.add_argument("--coloring_field",help="data field used for coloring node. Usual fields: diff, enrichment_subnet_score", default="diff")
    parser.add_argument("--SetScalarBarVisibility", help="show scalar bar or not", type=bool, default=False)
    parser.add_argument("--min_val",help="min value of legend bar", type=float, default=-1.0)
    parser.add_argument("--max_val",help="max value of legend bar", type=float, default=1.0)
    parser.add_argument("--ImageResolution_w",help="Image resolution, width", type=int, default=13990)
    parser.add_argument("--ImageResolution_h",help="Image resolution, height", type=int, default=8380)
    parser.add_argument("--opacity", help="node opacity", type=float, default=1.0)
    parser.add_argument("--delete_input_file", help="Delete vtu input file to save hard-disk mermory", type=bool, default=False)

    args = parser.parse_args()
    export_network_image(infile=args.infile,
                     outfile=args.outfile,
                     ImageResolution=[args.ImageResolution_w, args.ImageResolution_h],
                     coloring_field=args.coloring_field,
                     SetScalarBarVisibility=args.SetScalarBarVisibility,
                     min_val=args.min_val,
                     max_val=args.max_val,
                     opacity=args.opacity,
                    delete_input_file=args.delete_input_file)
