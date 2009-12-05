# name:     tidy-post.sed
# synopsis: script for postprocessing the output of
#           tidy(1) to get a valid XHTML document
# author:   Dr. Christoph L. Spiel


# The public id is in tidy.cfg.
1,9s|""|"http://www.w3.org/TR/MathML2/dtd/xhtml-math11-f.dtd"|

# Trim attibutes that are not allowed in XHTML.
s|<html\([^>]*\) lang="[^"]*"\([^>]*\)|<html\1\2|
s|<ol[^>]*>|<ol>|
s|<t\([dh]\)\([^>]*\)width="[^"]*"\([^>]*\)|<t\1\2\3|
s|<ul\([^>]*\)compact="[^"]*"\([^>]*\)|<ul\1\2|

# Undo our dirty work that tricks tidy(1) into correctly
# formatting inline math tags.
s|<mathinline|<math|g
s|</mathinline|</math|g
s|<mrowinline|<mrow|g
s|</mrowinline|</mrow|g
s|<msubinline|<msub|g
s|</msubinline|</msub|g
s|<msupinline|<msup|g
s|</msupinline|</msup|g
