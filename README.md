# LG-tangle-computer
**Author:** Matthew Harper

This repository provides standalone code which outputs the image of a given tangle under the Reshetikhin--Turaev functor for LG.

## Repository Description

The repository includes a Jupyter notebook `compute_LG.ipynb` and Python module `compute_LG_tangle.py`. The module constructs the Kac modules for quantized $$\mathfrak{sl}(2|1)$$ and the elementary tangle maps associated to the module $$V(0,\alpha)$$, then to any elementary tangle word it computes the associated RT morphism. The notebook provides an interface for running the module. One must provide a name for the tangle, a tangle word, and the width of the bottom of the tangle. When run it will display the output matrix, save a `name.pkl` for that matrix and a `name_computation.txt` which recalls the name, word, and bottom width. It also gives the elapsed time to compute and, for $$(1,1)$$-tangles, the scalar value of the matrix.
   
## Script

Only when run the very first time, build the elementary tangle maps by uncommenting the `%run compute_LG_tangle.py`. The convention for a tangle word is a comma delimited string with entries of the form `(number)+(elementary tangle)`. The `number` indicates the number of strands to the left of the elementary tangle. The elementary tangles are:

  - `er`, `el`: right and left evaluation
  - `cr`, `cl`: right and left coevaluation
  - `u`,`r`,`l`,`d`: positive crossings whose arrows point up, right, left, or down
  - `ui`,`ri`,`li`,`di`: negative crossings whose arrows point up, right, left, or down

The Hopf Link, for example, has tangle word: `"0cl, 1cl, 0u, 0u, 1er, 0er"`. The code will return the $$1\times 1$$ matrix $$0$$ since $$\mathrm{LG}$$ vanishes on closed diagrams.

## Warning

Since this computes matrices directly, one should try to use a tangle diagram with as small an index (width) as possible. A tangle diagram with maximum width $$d$$ will require a matrix with at least one dimension of size $$4^d$$. In my experience, one can expect to compute diagrams with $$d\leq 6$$ and in some cases $$d=7$$.

## Citation
If you use this code in your research, please cite this repository.

```bibtex
@misc{LG-tangle-computer,
  title = {LG-tangle-computer},
  author = {Harper, Matthew},
  year = {2026},
  note = {\url{https://github.com/mrhmath/LG-tangle-computer}},
}
```
## Acknowledgments
MH thanks Micah Chrisman, Ben-Michael Kohli, and Mathieu Kohli for helpful discussions. MH was partially supported through the NSF-RTG grant \#DMS-2135960.
