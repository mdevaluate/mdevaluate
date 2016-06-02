
Contribution guide
==================

This document aims to lay out the basics of contributing code to the ``mdevaluate`` package.
The code is managed through a git repository, hence this guides gives basic information on the usage of `git <https://git-scm.com>`_.
Int this document the prefix ``$`` indicates commands which should be ran on a shell.
For a brief 15 min interactive tutorial visit `try.github.org <https://try.gitbhub.org>`_.


Let's start with a short introduction to the terminology.
Python code is organized in *packages* and *modules*:

Modules:
  Any python file (e.g. ``test.py``) is called a module. A module can be imported (``import test``) an then used
  in other python code if it is in the python path, for example the working directory.
  In principle, importing a package means simply executing the code inside the file.
  All defintions, like variables or functions, are then available under the modules name.

Packages:
  Python modules can be grouped into packages. A python package is simply a folder,
  which contains at least one mandatory file ``__init__.py``. This file is the entry
  point into the module that is imported if the package is imported.
  All modules in the folder are treated as submodules, which can be accessed via
  a dot syntax, e.g. ``import package.test``. Packges can also contain subpackages.

A more `detailed explanation <https://docs.python.org/3/tutorial/modules.html>`_ can be found in the official python documentation.

Organization of the code
++++++++++++++++++++++++



The code for the evaluation software is organized in two python packges:

- ``pygmx``: This package provides a python wrapper for the Gromacs library and
  thereby functionality to read many of the file formats used within Gromacs.
- ``mdevaluate``: This package provides functionality for evaluation of molecular
  dynamics simulations. It uses the ``pygmx`` package to read files, but is
  (in theory) not  limited to Gromacs data.

Submodules
----------

Below the content of the submodules of the package is described.

atoms.py
........

Defintion of the ``Atom`` class and related stuff for atom selection and information.

autosave.py
...........

Experimental functionality for automatic saving and loading of evaluated data,
like correlation functions. For each function call a checksum is calculated
from the input, which only changes if the output of the function changes.

coordinates.py
..............

Definition of the ``Coordinates`` class and ``CoordinatesMap`` for coordinates
transformations and related functions.

correlation.py
..............

Functionality to evaluate correlation functions.

distribution.py
...............

Functionality to evaluate distribution functions.

utils.py
........

A collection of utility functions.

Set up a development environment
++++++++++++++++++++++++++++++++

The code is hosted in the groups git server https://chaos3.fkp.physik.tu-darmstadt.de,
use your credentials of the institutes computers for login.
Authentication has to be done through SSH-Keys, therefore first add the SSH Key of
your development machine to your account.

1. Copy the public key from the file ``~/.ssh/id_rsa.pub``
2. Add it to your account at https://chaos3.fkp.physik.tu-darmstadt.de/settings/panel/ssh/
3. Clone the repository

::

  $ git clone ssh://vcs@chaos3.fkp.physik.tu-darmstadt.de/diffusion/MDE/mdevaluate.git

Organization of the repository
------------------------------

The repository is organized through git branches.
At the moment there exist two branches in the remote repository: *master* and *dev*.


Adding code to the repository
+++++++++++++++++++++++++++++

All changes to the code are done in your local clone of the repository.
If a feature is complete, or at least works, the code can be pushed to the remote,
to make it accessible for others.

A standard workflow to submit new code is the following

1. Create a new **branch** in your local repository
2. **Commit** the changes to your local repository
3. **Merge** the branch into master
4. **Push** the changes to the remote reposiory

Pulling updates from remote
---------------------------

Before working with the code, the latest updates should be pulled for the master branch::

  $ git checkout master
  $ git pull

Create a new branch
-------------------

Before changing any code, create a new branch in your local repository.
This helps to keep an overview of all the changes and simplifies merging.
To create a new branch locally enter the following commands::

  $ git checkout master
  $ git branch my-feature
  $ git checkout my-feature

First switch to the master branch to make sure the new branch is based on it.
Then create the new branch, called `my-fetaure` and switch to it.
Now you can start making changes in the code.

Commiting changes
-----------------

A bundle of changes in the code is called a *commit*.
These changes can happen in different files and should be associated with each other.
Let's assume, two files have been changed (``atoms.py`` and ``utils.py``).
The command::

  $ git diff atoms.py

will show you all changes that were made in the file since the latest commit.
Before commiting changes have to be *staged*, which is done by::

  $ git add atoms.py utils.py

This my be repeated as often as necessary.
When all changes for a commit are staged, it can actually be created::

  $ git commit

This will open up an editor where a commit message has to be entered.
After writing the commit message, save & close the file, which will create the commit.

Merging into master
-------------------

When all changes are made and the new feature should be made public, first the branch has to be merged into master.
Most of the time, the master branch will have been updated, therfore first pull any updates::

  $ git checkout master
  $ git pull

When the master branch is up to date, it can be merged into the feature branch::

  $ git merge my-feature

If no conflicting changes were made, merging works automatically.
If for example the same line was modified in a commit in master and your commits, a merge conflict will occur.
Git tells you which files have conflicts and asks you to resolve these.
The respective lines will be marked with conflict-resolution markers in the files.
The most basic way of resolving a conflivt is by editing these files and choosing the appropiate version of the code.
See the `git documentation <https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging#Basic-Merge-Conflicts>`_ for an explanation.
After resolving the conflict, the files need to be staged and the merge has to be committed::

  $ git add utils.py
  $ git commit

The commit message will be generated automatically, indicating the merge.

Push to remote
--------------

After merging the changes can be pushed to the remote::

  $ git push

The new code is now available in the remote.
