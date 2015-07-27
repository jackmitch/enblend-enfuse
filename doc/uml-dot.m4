m4_dnl  Definitions for a common structure and appearance of all flow charts


m4_dnl  Attributes of the `Start' and `Finish' symbols.

m4_define(`uml_terminal_color', `black')
m4_define(`uml_terminal_size', `0.2')

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


m4_dnl  Shape of Activity boxes

m4_define(`uml_activity', `shape = box, style = rounded')


m4_dnl  Shape of Notes and their associated edges

m4_define(`uml_note_color', `yellow')

m4_define(`uml_note', `fillcolor = uml_note_color, shape = note, style = filled')
m4_define(`uml_note_edge', `dir = none, style = dashed')


m4_dnl  Shape and contents of conditionals

m4_define(`uml_terminal_height', `0.3')
m4_define(`uml_terminal_width', `0.6')

m4_define(`uml_branch',
          `fixedsize = shape,
           height = uml_terminal_height,
           label = "",
           shape = diamond,
           style = solid,
           width = uml_terminal_width')
m4_define(`uml_merge', `uml_branch')
