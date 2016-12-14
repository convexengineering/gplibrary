Testing of external repositories
================================

Directions
----------

  2. In your repository's top level, create the file `TESTS`.
     Each line should be the local path to a `.py` file containing a `test()` function, like this:
     
         mission.py
         aircraft/wing/wing_spar.py

  3. (optional) create a `TESTCONFIG` file in the same folder. The following options are supported:
      `pip install` specifies packages your repository requires,
      `gpkit-models` the branch of gpkit-models you're working on,
      and `skipsolvers` the solvers that cannot solve your code.
      Format the file like this:
  
        pip install : pandas
        gpkit-models branch : 1682
        skipsolvers : cvxopt

  4. Test your repository locally by running `python -c "from gpkit.tests.from_paths import run; run()"`
    (consider saving this line as `runtests.sh` for easy access)
    
  1. Submit a PR to add the name of your repository to the `EXTERNALTESTS` file in `gpkit-models`
