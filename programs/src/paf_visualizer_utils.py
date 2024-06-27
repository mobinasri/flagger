from block_utils import *
from alignment_utils import *

def readAnnotationsWithColors(bedPath, nameColIndex, colorColIndex, isColorRGB, isEightBit, patternForName):
  """

  bedPath: Path to the annotation bed file
  nameColIndex: 0-based index of the column in bed file that contains the name of the annotation
  colorColIndex: 0-based index of the column in bed file that contains the related color
  isColorRGB: True if the color is given as a three numbers sperated by comma representing the RGB value of the color
              it can be given as either 8 bit values like "125,200,80" or ratios like "0.1,0.5,0.8"
  isEightBit: True if the RGB color given as either 8 bit values like "125,200,80"
  patternForName: A pattern for extracting the name of each track. It is useful when there are some other information
                  in the name column that we don't want save in the output list. For example in a censat annotation
                  we might have "dhor(S01CMH1d)" but we only want the "dhor" part.
  """
  annotBlockListPerContig = BlockList.parseBed(bedPath, saveAllOtherColumns=True)
  for contig, blockList in annotBlockListPerContig.items():
    for i in range(len(blockList.blocks)):
      block = blockList.blocks[i]
      color = block[2][colorColIndex - 3]
      # if color is given as something like "125,200,80"
      if isColorRGB and isEightBit:
        color = [int(c)/ 255 for c in color.split(",")]
      # if color is given as something like "0.1,0.5,0.8"
      if isColorRGB and not isEightBit:
        color = [float(c) for c in color.split(",")]

      name = block[2][nameColIndex - 3]
      # extract a name based on the given pattern
      if patternForName != '' and patternForName != None:
        matchedName='NA'
        m = re.search(patternForName, name)
        if m:
          matchedName = m.group(0)
      else:  # if no pattern is given use the name itself
        matchedName = name
      # save the third entry as dict of name and color
      blockList.blocks[i] = (block[0], block[1], {"name": matchedName, "color": color})
  return annotBlockListPerContig



def getIntervalEndPoints(alignment):
  rStart = alignment.chromStart
  rEnd = alignment.chromEnd
  rLen = alignment.chromLength
  qStart = alignment.contigStart
  qEnd = alignment.contigEnd
  qLen = alignment.contigLength

  x = rStart
  y = qEnd if alignment.orientation == '-' else qStart
  xPoints = [x]
  yPoints = [y]
  for op, opSize in alignment.cigarList:
    if op == 'X' or op == '=':
      if alignment.orientation == '+':
        x += opSize
        y += opSize
      else:
        x += opSize
        y -= opSize
    if op == 'I':
      if alignment.orientation == '+':
        y += opSize
      else:
        y -= opSize
    if op == 'D':
      if alignment.orientation == '+':
        x += opSize
      else:
        x += opSize

    xPoints.append(x)
    yPoints.append(y)
  return np.array(xPoints), np.array(yPoints)


plotAttributes = {'dotplot': {'lineWidth' : 1,
                              'xCoorsInMb': True,
                              'yCoorsInMb': True,
                              'xLim': (10e6, 20e6),
                              'yLim': (10e6, 20e6),
                              'plotBreakPoints': True},
                  'breakPoints': {'ratioOfPanel': 0.05,
                                  'color': 'red',
                                  'fontSize': 8},
                  'annotation': {'ratioOfPanel': 0.05,
                                 'fontSize': 8}}

def plotAlignments(alignments, plotAttributes, refAnnotBlockList=[], queryAnnotBlockList=[], refAnnotNames = [], queryAnnotNames=[]):

  panels = []
  fig = plt.figure(figsize=(8,8))
  axMainPanel = plt.axes([0.1,0.1, 0.6,0.6])
  panels.append(axMainPanel)


  # Set scales if x and y coordinates should be in Mb
  xScale = 1.0
  yScale = 1.0
  if plotAttributes['dotplot']['xCoorsInMb'] :
    xScale = 1.0e6
  if plotAttributes['dotplot']['yCoorsInMb'] :
    yScale = 1.0e6

  # Each Panel contains two x/y limits
  # The inner limits; (xLimMinDotPlot, xLimMaxDotPlot) and (yLimMinDotPlot, yLimMaxDotPlot)
  # are for the dotplot of the alignments parsed from the paf file
  # The outer limits; (xLimMinPanel, xLimMaxPanel) and (yLimMinPanel, yLimMaxPanel)
  # are for the whole panel, which contains the inner dotplot and
  # also break points and annotation tracks if given

  xLimMaxDotPlot = plotAttributes['dotplot']['xLim'][1] / xScale
  xLimMinDotPlot = plotAttributes['dotplot']['xLim'][0] / xScale
  xSpanDotPlot = xLimMaxDotPlot - xLimMinDotPlot


  yLimMaxDotPlot = plotAttributes['dotplot']['yLim'][1] / yScale
  yLimMinDotPlot = plotAttributes['dotplot']['yLim'][0] / yScale
  ySpanDotPlot = yLimMaxDotPlot - yLimMinDotPlot


  # The panel limits are initialized
  yLimMinPanel = yLimMinDotPlot
  yLimMaxPanel = yLimMaxDotPlot

  xLimMinPanel = xLimMinDotPlot
  xLimMaxPanel = xLimMaxDotPlot

  # The panel limits are extended if any annotation is given
  # The extension are based on the number of given annotation
  # and the amount space each annotation should take

  # Extending panel limits on the y axis if ref annotation is given
  if refAnnotBlockList != []:
    refAnnotHeight = plotAttributes['annotation']['ratioOfPanel'] * ySpanDotPlot
    yLimMinPanel -= refAnnotHeight * len(refAnnotBlockList)
  else:
    refAnnotHeight = 0

  # Extending panel limits on the x axis if query annotation is given
  if queryAnnotBlockList != []:
    queryAnnotWidth = plotAttributes['annotation']['ratioOfPanel'] * xSpanDotPlot
    xLimMinPanel -= queryAnnotWidth * len(queryAnnotBlockList)
  else:
    queryAnnotWidth = 0

  # Extend the panel limits more if the alignment break points should be plotted
  # alignment break points mean start/end of each alignment
  if plotAttributes['dotplot']['plotBreakPoints']:
    breakPointHeight = plotAttributes['breakPoints']['ratioOfPanel'] * ySpanDotPlot
    breakPointWidth = plotAttributes['breakPoints']['ratioOfPanel'] * xSpanDotPlot
    xLimMinPanel -= breakPointWidth
    yLimMinPanel -= breakPointHeight
  else:
    breakPointWidth = 0
    breakPointHeight = 0
  #
  # -------------------------------
  # |                             |
  # |    ---------------------    |
  # |    |                    |   |
  # |    |                    |   |
  # |    |                    |   |
  # |    |   DotPlot Panel    |   |
  # |    |                    |   |  <-- Whole Panel
  # |    |                    |   |
  # |    |                    |   |
  # |    |                    |   |
  # |     --------------------    |
  # |                             |
  # ------------------------------


  # xMargin and yMargin are defined here
  # xMargin (for x axis) or yMargin (for y axis) is the size of the gap between the first (outer)
  # annotation track and the border of the whole panel. It is also the size of the gap between
  # the break point track and the border of the dotplot panel.
  # Between annotation tracks there is a gap of xMagin/2 (for x axis) or yMargin/2 (for y axis)
  # Similarly between the last annotatoin track and the break point track there is a gap of
  # xMagin/2 (for x axis) or yMargin/2 (for y axis)
  yMargin = 0
  xMargin = 0
  if breakPointHeight > 0 or refAnnotHeight > 0:
    yMargin = max(breakPointHeight, refAnnotHeight)
  if breakPointWidth > 0 or queryAnnotWidth > 0:
    xMargin = max(breakPointWidth, queryAnnotWidth)

  # 2 = 1 (between annotation and the whole panle) + 1 (between break point and the dotplot panel)
  xMarginFactor = (2 + 0.5 * len(queryAnnotBlockList))
  yMarginFactor = (2 + 0.5 * len(refAnnotBlockList))
  xLimMinPanel -= xMargin * xMarginFactor
  yLimMinPanel -= yMargin * yMarginFactor

  xWholeLeftMargin = xMargin * xMarginFactor + queryAnnotWidth * len(queryAnnotBlockList) + breakPointWidth
  yWholeBottomMargin = yMargin * yMarginFactor + refAnnotHeight * len(refAnnotBlockList) + breakPointHeight

  xWholeRightMargin = xWholeLeftMargin
  yWholeTopMargin = yWholeBottomMargin

  xLimMaxPanel += xWholeRightMargin
  yLimMaxPanel += yWholeTopMargin



  xPointMax = 0
  yPointMax = 0
  for alignment in alignments:
    if (alignment.chromEnd / xScale) < xLimMinDotPlot and (alignment.contigEnd / yScale) < yLimMinDotPlot:
      continue
    if (alignment.chromStart / xScale) > xLimMaxDotPlot and (alignment.contigStart / yScale) > yLimMaxDotPlot:
      continue
    #if alignment.contigEnd - alignment.contigStart < 50e3:
    #  continue
    xPoints, yPoints = getIntervals(alignment)

    xPoints = xPoints / xScale
    yPoints = yPoints / yScale

    # keep only the points in the requested x/y limits
    keepIndices = (xLimMinDotPlot < xPoints) &  (xPoints < xLimMaxDotPlot) & (yLimMinDotPlot < yPoints) & (yPoints < yLimMaxDotPlot)
    # skip it if there is no part of the alignment intersecting with the limits
    if(sum(keepIndices) == 0):
      continue
    xPoints = xPoints[keepIndices]
    yPoints = yPoints[keepIndices]

    axMainPanel.plot(xPoints, yPoints, linewidth=plotAttributes['dotplot']['lineWidth'], zorder=1)
    xPointMax = max(xPointMax, max(xPoints))
    yPointMax = max(yPointMax, max(yPoints))


    #####################
    # Plot break points #
    #####################

    # plot if the related attribute is set
    if plotAttributes["dotplot"]["plotBreakPoints"]:

      # Reference start point on x axis
      if xLimMinDotPlot < (alignment.chromStart / xScale) < xLimMaxDotPlot:
        axMainPanel.plot([alignment.chromStart / xScale,
                alignment.chromStart / xScale],
                [yLimMinDotPlot - yMargin - breakPointHeight,
                yLimMinDotPlot - yMargin],
                color=plotAttributes['breakPoints']['color'],
                alpha=plotAttributes['breakPoints']['alpha'],
                linestyle="solid",
                linewidth=1,
                zorder=3)

      # Reference end point on x axis
      if xLimMinDotPlot < (alignment.chromEnd / xScale) < xLimMaxDotPlot:
        axMainPanel.plot([alignment.chromEnd / xScale,
                 alignment.chromEnd / xScale],
                [yLimMinDotPlot - yMargin - breakPointHeight,
                 yLimMinDotPlot - yMargin],
                color = plotAttributes['breakPoints']['color'],
                alpha=plotAttributes['breakPoints']['alpha'],
                linestyle="solid",
                linewidth=1,
                zorder = 3)

      # Query start point on y axis
      if yLimMinDotPlot < (alignment.contigStart / yScale) < yLimMaxDotPlot:
        axMainPanel.plot([xLimMinDotPlot - xMargin - breakPointWidth,
                xLimMinDotPlot - xMargin],
                [alignment.contigStart / yScale,
                alignment.contigStart / yScale],
                color=plotAttributes['breakPoints']['color'],
                alpha=plotAttributes['breakPoints']['alpha'],
                linestyle="solid",
                linewidth=1,
                zorder=3)

      # Query end point on y axis
      if yLimMinDotPlot < (alignment.contigEnd / yScale) < yLimMaxDotPlot:
        axMainPanel.plot([xLimMinDotPlot - xMargin - breakPointWidth,
                 xLimMinDotPlot - xMargin],
                [alignment.contigEnd / yScale,
                 alignment.contigEnd / yScale],
                color=plotAttributes['breakPoints']['color'],
                alpha=plotAttributes['breakPoints']['alpha'],
                linestyle="solid",
                linewidth=1,
                zorder=3)

  # set the limits of the whole panel
  axMainPanel.set_xlim((xLimMinPanel, xLimMaxPanel))
  axMainPanel.set_ylim((yLimMinPanel, yLimMaxPanel))


  #################################################
  # Write "Brk" at the ends of break point tracks #
  #################################################

  if plotAttributes["dotplot"]["plotBreakPoints"]:
    axMainPanel.text(s="Brk",
            x = xLimMaxDotPlot + xWholeRightMargin * 0.1,
            y = yLimMinDotPlot - yMargin - breakPointHeight / 2,
            horizontalalignment = 'left',
            verticalalignment = 'center',
            fontsize = plotAttributes['breakPoints']['fontSize'],
            zorder=3)
    axMainPanel.text(s="Brk",
            x = xLimMinDotPlot - xMargin - breakPointWidth/2,
            y = yLimMaxDotPlot + yWholeTopMargin * 0.1,
            horizontalalignment = 'center',
            verticalalignment = 'bottom',
            fontsize = plotAttributes['breakPoints']['fontSize'],
            rotation = 90,
            zorder=3)

  ###########################################################
  # Write annotation names at the ends of annotation tracks #
  ###########################################################

  # Reference annotation along the x axis
  for i in range(len(refAnnotNames)):
    axMainPanel.text(s = refAnnotNames[i],
            x = xLimMaxDotPlot + xWholeRightMargin * 0.1,
            y = yLimMinPanel + yMargin * (1 + 0.5 * i) + refAnnotHeight * (i + 0.5),
            horizontalalignment = 'left',
            fontsize = plotAttributes['annotation']['fontSize'],
            verticalalignment = 'center',
            zorder=3)

  # Query annotation along the y axis
  for i in range(len(queryAnnotNames)):
    axMainPanel.text(s = queryAnnotNames[i],
            x = xLimMinPanel + xMargin * (1 + 0.5 * i) + queryAnnotWidth * (i + 0.5),
            y = yLimMaxDotPlot + yWholeTopMargin * 0.1,
            horizontalalignment = 'center',
            verticalalignment = 'bottom',
            fontsize = plotAttributes['annotation']['fontSize'],
            rotation = 90,
            zorder=3)


  existingAnnotations = {}

  ##########################
  # Plot annotation tracks #
  ##########################

  # annotation tracks for the query sequence on y axis
  annotRects = []
  if queryAnnotBlockList != []:
    for i in range(len(queryAnnotBlockList)):
      existingAnnotations[queryAnnotNames[i]] = {}
      for block in queryAnnotBlockList[i].blocks:
        # skip the annotation block if it does not have any overlap
        # with the  x/y limits
        if (block[1] / yScale) < yLimMinDotPlot or (block[0] / yScale) > yLimMaxDotPlot:
          continue
        # plot only part of the annotation overlapping the x/y limits
        s = max(block[0] / yScale, yLimMinDotPlot)
        e = min(block[1] / yScale, yLimMaxDotPlot)
        annotRect = Rectangle((xLimMinPanel + xMargin * (1 + 0.5 * i) + queryAnnotWidth * i, s),
                               width = queryAnnotWidth,
                               height = e-s,
                               color=block[2]["color"])
        existingAnnotations[queryAnnotNames[i]][block[2]["name"]] = block[2]["color"]
        annotRects.append(annotRect)

  # annotation tracks for the ref sequence on x axis
  if refAnnotBlockList != []:
    for i in range(len(refAnnotBlockList)):
      existingAnnotations[refAnnotNames[i]] = {}
      for block in refAnnotBlockList[i].blocks:
        # skip the annotation block if it does not have any overlap
        # with the  x/y limits
        if (block[1] / xScale) < xLimMinDotPlot or (block[0] / xScale) > xLimMaxDotPlot:
          continue
        # plot only part of the annotation overlapping the x/y limits
        s = max(block[0] / xScale, xLimMinDotPlot)
        e = min(block[1] / xScale, xLimMaxDotPlot)
        annotRect = Rectangle((s, yLimMinPanel + yMargin * (1 + 0.5 * i) + refAnnotHeight * i),
                               width = e-s,
                               height = refAnnotHeight,
                               color=block[2]["color"])
        existingAnnotations[refAnnotNames[i]][block[2]["name"]] = block[2]["color"]
        annotRects.append(annotRect)

  axMainPanel.add_collection(PatchCollection(annotRects, match_original=True, zorder=3))

  ###############################
  # Plot lenged for annotations #
  ###############################

  widthRatioForLegendPanel = 0.2
  heightRatioForLegendPanel = 0.2
  gapRatioForLegendPanel = 0.05
  mainPanelTopRatio = 0.7
  mainPanelRightRatio = 0.7
  j=1
  for annotName, annotDict in existingAnnotations.items():
    axLegend = plt.axes([mainPanelRightRatio * 1.05,
                          mainPanelTopRatio - j * heightRatioForLegendPanel - (j-1) * gapRatioForLegendPanel,
                          widthRatioForLegendPanel,
                          heightRatioForLegendPanel])
    axLegend.set_yticks([])
    axLegend.set_yticklabels([])
    axLegend.set_xticks([])
    axLegend.set_xticklabels([])
    panels.append(axLegend)
    axLegend.text(s = annotName,
                  x = 0.5,
                  y = 0.95,
                  horizontalalignment='center',
                  verticalalignment='top',
                  fontsize = 8)

    x = 0.1
    y = 0.1
    annotRectWidth = 0.3
    annotRectHeight = 0.7 / len(annotDict) * 0.7
    annotRectGap = 0.7 / len(annotDict) * 0.3
    for name, color in annotDict.items():
      annotRect = Rectangle((x, y),
                             width = annotRectWidth,
                             height = annotRectHeight,
                             color=color)
      axLegend.add_patch(annotRect)
      axLegend.text(s=name,
                    x = 0.75,
                    y = y + annotRectHeight / 2,
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=8)
      y += annotRectHeight + annotRectGap
    j += 1



  # plot a border around the dotplot panel
  # slightly larger than the panel (by 0.005)
  borderRect = Rectangle((xLimMinDotPlot - xSpanDotPlot * 0.005, yLimMinDotPlot - ySpanDotPlot* 0.005),
                         width = xSpanDotPlot * 1.01 ,
                         height = ySpanDotPlot * 1.01,
                         facecolor='none',
                         alpha=1,
                         edgecolor='gray',
                         linewidth=1,
                         zorder = 4)
  axMainPanel.add_patch(borderRect)
  return fig, panels
