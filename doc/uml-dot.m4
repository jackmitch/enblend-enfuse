m4_dnl  This file is part of Enblend.
m4_dnl  Licence details can be found in the file COPYING.


m4_dnl  Definitions for a common structure and appearance of all
m4_dnl  flow charts.
m4_dnl
m4_dnl  Despite the name of the file and the prefix of all macros
m4_dnl  (`uml_'), we treat only `Activity Diagrams' here for we do not
m4_dnl  need more for the Enblend/Enfuse documentation.


m4_define(`uml_font', `Helvetica')


m4_dnl  Graph Attributes

m4_define(`uml_graph_font', `uml_font')
m4_define(`uml_graph_font_size', `10')

m4_define(`uml_activity_graph',
          `fontname = uml_graph_font,
           fontsize = uml_graph_font_size,
           forcelabels = true,
           splines = ortho')

m4_dnl _Sometimes we want to fit a larger graph to a printed page.
m4_dnl  This macro helps while leaving other formats than EPS anone.
m4_define(`uml_compressed_layout',
          m4_ifdef(`dot_output_type',
                   m4_ifelse(dot_output_type, `eps',
                             `ranksep = 0.375, ratio = compress')))


m4_dnl  Attributes of the `Start' and `Finish' terminal symbols

m4_define(`uml_terminal_color', `black')
m4_define(`uml_terminal_size', `0.125')

m4_define(`uml_start',
          `fillcolor = uml_terminal_color,
           fixedsize = shape,
           label = "",
           shape = circle,
           style = filled,
           width = uml_terminal_size')

m4_define(`uml_finish',
          `fillcolor = uml_terminal_color,
           fixedsize = shape,
           label = "",
           shape = doublecircle,
           style = filled,
           width = uml_terminal_size')


m4_dnl  Shape of Activity boxes and the edges that connect them

m4_define(`uml_activity_font', `uml_font')
m4_define(`uml_activity_font_size', `9')
m4_define(`uml_activity_penwidth', `0.5')
m4_define(`uml_arrow_size', `0.667')
m4_define(`uml_edge_font_size', `8')

m4_define(`uml_activity',
          `fontname = uml_activity_font,
           fontsize = uml_activity_font_size,
           penwidth = uml_activity_penwidth,
           shape = box,
           style = rounded')
m4_define(`uml_edge',
          `arrowsize = uml_arrow_size,
           fontname = uml_activity_font,
           fontsize = uml_edge_font_size,
           penwidth = uml_activity_penwidth')


m4_dnl  Shape of Notes and the edges that connect Activity boxes and Notes

m4_define(`uml_note_color', `"0.167,0.4,1.0"') # RGB: 0xffff99, "pale yellow"
m4_define(`uml_note_font', `uml_font')
m4_define(`uml_note_font_size', `9')
m4_define(`uml_note_penwidth', `0.35')

m4_define(`uml_note',
          `fillcolor = uml_note_color,
           fontname = uml_note_font,
           fontsize = uml_note_font_size,
           penwidth = uml_note_penwidth,
           shape = note,
           style = filled')
m4_define(`uml_note_edge',
          `dir = none,
           fontname = uml_note_font,
           fontsize = uml_edge_font_size,
           penwidth = uml_note_penwidth,
           style = dashed')


m4_dnl  Shape and contents of Branch/Merge conditionals

m4_define(`uml_conditional_height', `0.25')
m4_define(`uml_conditional_width', `0.5')

m4_define(`uml_branch',
          `fixedsize = shape,
           height = uml_conditional_height,
           label = "",
           shape = diamond,
           style = solid,
           width = uml_conditional_width')
m4_define(`uml_merge', `uml_branch')


m4_dnl  Shape and contents of Join/Fork parallelism symbols

m4_define(`uml_parallelism_height', `0.02')
m4_define(`uml_parallelism_width', `2')

m4_define(`uml_fork',
          `color = black,
           fixedsize = shape,
           height = uml_parallelism_height,
           label = "",
           shape = box,
           style = filled,
           width = uml_parallelism_width')
m4_define(`uml_join', `uml_fork')
