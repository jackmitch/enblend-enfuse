m4_dnl  This file is part of Enblend.
m4_dnl  Licence details can be found in the file COPYING.


m4_dnl  Definitions for a common structure and appearance of all
m4_dnl  flow charts.
m4_dnl
m4_dnl  Despite the name of the file and the prefix of all macros
m4_dnl  (`uml_'), We treat only `Activity Diagrams' here.  We do not
m4_dnl  need more for the Enblend/Enfuse documentation.


m4_dnl  Graph Attributes

m4_define(`uml_activity_graph', `forcelabels = true')


m4_dnl  Attributes of the `Start' and `Finish' symbols.

m4_define(`uml_terminal_color', `black')
m4_define(`uml_terminal_size', `0.167')

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


m4_dnl  Shape of Activity boxes and their associated edges

m4_define(`uml_activity_penwidth', `0.5')
m4_define(`uml_arrow_size', `0.667')

m4_define(`uml_activity',
          `penwidth = uml_activity_penwidth,
           shape = box,
           style = rounded')
m4_define(`uml_edge',
          `arrowsize = uml_arrow_size,
           penwidth = uml_activity_penwidth')


m4_dnl  Shape of Notes and their associated edges

m4_define(`uml_note_color', `"0.167,0.4,1.0"') # RGB: 0xffff99
m4_define(`uml_note_penwidth', `0.35')

m4_define(`uml_note',
          `fillcolor = uml_note_color,
           penwidth = uml_note_penwidth,
           shape = note,
           style = filled')
m4_define(`uml_note_edge',
          `dir = none,
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
