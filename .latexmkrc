@default_excluded_files = ( 'mpcPackageJulia_glos.tex' );
$pdf_mode = 4;
$lualatex = "lualatex -synctex=1 -interaction=nonstopmode -file-line-error -shell-escape %O %S";