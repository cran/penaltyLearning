(TeX-add-style-hook
 "Definition"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fullpage" "cm")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "fullpage"
    "verbatim"
    "hyperref"
    "graphicx"
    "natbib"
    "amsmath"
    "amssymb")
   (TeX-add-symbols
    "argmin"
    "Diag"
    "TPR"
    "FPR"
    "FN"
    "FP"
    "argmax"
    "maximize"
    "minimize"
    "RR")
   (LaTeX-add-labels
    "eq:seg-model"
    "eq:modelSelection"))
 :latex)

