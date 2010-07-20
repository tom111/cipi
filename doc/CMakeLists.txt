ADD_CUSTOM_COMMAND(
  OUTPUT    ${CMAKE_CURRENT_BINARY_DIR}/manual.aux
  DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/manual.tex
  COMMAND   pdflatex
  ARGS      -interaction=batchmode ${CMAKE_CURRENT_BINARY_DIR}/manual
  COMMENT   "Latex (first pass)"
  )

ADD_CUSTOM_COMMAND(
  OUTPUT    ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf
  DEPENDS   ${CMAKE_CURRENT_BINARY_DIR}/manual.aux
  COMMAND   pdflatex
  ARGS      -interaction=batchmode ${CMAKE_CURRENT_BINARY_DIR}/manual
  COMMENT   "Latex (second pass)"
  )

# Eventually trigger the whole process
ADD_CUSTOM_TARGET(LaTeXDocument ALL echo
   DEPENDS   ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf
   )

INSTALL(FILES manual.pdf DESTINATION share/doc/cipi-${CIPI_VERSION})