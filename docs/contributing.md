# Contributors Guide

The adjTomo software suite is community driven project and everyone is welcome 
to contribute to any of the packages. The source code is hosted on 
[GitHub](https://github.com/adjtomo/pyatoa).

*This guide is modified from 
[PyGMTs](https://github.com/GenericMappingTools/pygmt) Contributors guide*

## Ways to Contribute

### Ways to Contribute Documentation and/or Code

* Tackle any issue that you wish! Issues labeled **"good first issue"** 
  indicate that they are beginner friendly, meaning that they don't require 
  extensive knowledge of the project.
* Make a tutorial or gallery example of how to do something.
* Improve the API documentation.
* Contribute code! This can be code that you already have and it doesn't need to  
  be perfect! We will help you clean things up, test it, etc.

### Ways to Contribute Feedback

* Provide feedback about how we can improve the project or about your particular use
  case. Open an issue, feature request or bug fix.
* Help solve issues, or give a "thumbs up" on issues that others reported which are
  relevant to you.

## Providing Feedback

### Reporting a Bug

* Find the [*Issues*](https://github.com/adjtomo/pyatoa/issues) tab on the
top of the GitHub repository and click *New issue*.
* Click on *Get started* next to *Bug report*.
* After submitting your bug report, try to answer any follow up questions about the bug
  as best as you can.

### Submitting a Feature Request

* Find the [*Issues*](https://github.com/adjtomo/pyatoa/issues) tab on the
  top of the GitHub repository and click *New issue*.
* Click on *Get started* next to *Feature request*.
* After submitting your feature request, try to answer any follow up questions as best
  as you can.

### Submitting General Comments/Questions

There are several pages on the [Community Forum](https://github.com/adjtomo/discussions)
where you can submit general comments and/or questions:

## General Guidelines

### Resources for New Contributors

Please take a look at these resources to learn about Git and pull requests (don't
hesitate to ask questions.

* [How to Contribute to Open Source](https://opensource.guide/how-to-contribute/).
* [Git Workflow Tutorial](http://www.asmeurer.com/git-workflow/) by Aaron Meurer.
* [How to Contribute to an Open Source Project on GitHub](https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github).

### Pull Request Workflow

We follow the [git pull request workflow](http://www.asmeurer.com/git-workflow)
to make changes to our codebase. Every change made goes through a pull request, even
our own, so that our
[continuous integration](https://en.wikipedia.org/wiki/Continuous_integration)
services have a chance to check that the code is up to standards and passes all
our tests. This way, the *main* branch is always stable.

#### General Guidelines for Making a Pull Request (PR):

* What should be included in a PR
  - Have a quick look at the titles of all the existing issues first. If there
    is already an issue that matches your PR, leave a comment there to let us
    know what you plan to do. Otherwise, **open an issue** describing what you
    want to do.
  - Each pull request should consist of a **small** and logical collection of
    changes; larger changes should be broken down into smaller parts and
    integrated separately.
  - Bug fixes should be submitted in separate PRs.
* How to write and submit a PR
  - Use underscores for all Python (*.py) files as per
    [PEP8](https://www.python.org/dev/peps/pep-0008/), not hyphens. Directory
    names should also use underscores instead of hyphens.
  - Describe what your PR changes and *why* this is a good thing. Be as
    specific as you can. The PR description is how we keep track of the changes
    made to the project over time.
  - Do not commit changes to files that are irrelevant to your feature or
    bugfix (e.g.: `.gitignore`, IDE project files, etc).
  - Write descriptive commit messages. Chris Beams has written a
    [guide](https://chris.beams.io/posts/git-commit/) on how to write good
    commit messages.
* PR review
  - Be willing to accept criticism and work on improving your code; we don't
    want to break other users' code, so care must be taken not to introduce
    bugs.
  - Be aware that the pull request review process is not immediate, and is
    generally proportional to the size of the pull request.

#### General Process for Pull Request Review:

After you've submitted a pull request, you should expect to hear at least a
comment within a couple of days. We may suggest some changes, improvements or
alternative implementation details.

To increase the chances of getting your pull request accepted quickly, try to:

* Submit a friendly PR
  - Write a good and detailed description of what the PR does.
  - Write some documentation for your code (docstrings) and leave comments
    explaining the *reason* behind non-obvious things.
  - Write tests for the code you wrote/modified if needed.
  - Include an example of new features in the gallery or tutorials.
* Have a good coding style
  - Use readable code, as it is better than clever code (even with comments).
  - Follow the [PEP8](http://pep8.org) style guide for code and the
    [NumPy style guide](https://numpydoc.readthedocs.io/en/latest/format.html)
    for docstrings. Please refer to [Code style](contributing.md#code-style).


## Contributing Documentation

### Documentation Overview

The documentation is written in
[reStructuredText](https://docutils.sourceforge.io/rst.html) and built by
[Sphinx](http://www.sphinx-doc.org/) and hosted on 
[ReadTheDocs](https://readthedocs.org). When contributing documentation, please
follow the general guidelines in the [pull request workflow](contributing.md#pull-request-workflow)
section.

Additionally, some documentation pages are written as Jupyter notebooks and 
converted to reStructuredText files with a Python script.

There are two primary ways to edit the Pyatoa documentation:
- For simple documentation changes, you can easily
  [edit the documentation on GitHub](contributing.md#editing-the-documentation-on-github).
  This only requires you to have a GitHub account.
- For more complicated changes, you can
  [edit the documentation locally](contributing.md#editing-the-documentation-locally).
  In order to build the documentation locally, you first need to
  [set up your environment](contributing.md#setting-up-your-environment).

### Editing the Documentation on GitHub

If you're browsing the documentation and notice a typo or something that could be
improved, please consider letting us know by [creating an issue](contributing.md#reporting-a-bug) or
(even better) submitting a fix.

You can submit fixes to the documentation pages completely online without having to
download and install anything:

1. On each documentation page, there should be an "Edit on GitHub" link at the very
  top.
2. Click on that link to open the respective source file on GitHub for editing 
   online (you'll need a GitHub account).
3. Make your desired changes.
4. When you're done, scroll to the bottom of the page.
5. Fill out the two fields under "Commit changes": the first is a short title describing
  your fixes; the second is a more detailed description of the changes. Try to be as
  detailed as possible and describe *why* you changed something.
6. Choose "Create a new branch for this commit and start a pull request" and
  click on the "Propose changes" button to open a pull request.
8. We'll review your pull request, recommend changes if necessary, and then merge
  them in if everything is OK.
9. Done!

Alternatively, you can make the changes offline to the files in the `doc` folder or the
example scripts. See [editing the documentation locally](contributing.md#editing-the-documentation-locally)
for instructions.

### Editing the Documentation Locally

For more extensive changes, you can edit the documentation in your cloned repository
and build the documentation to preview changes before submitting a pull request.

First you'll have to set up a Conda environment with the correct pacakges for
building documentation.

```bash
cd docs/
conda env create --file environment.yml
conda activate pyatoa-docs
```

After making your changes, you can build the HTML files from sources using:

```bash
make html
```

This will build the HTML files in `docs/_build/html`.
Open `doc/_build/html/index.html` in your browser to view the pages. Follow the
[pull request workflow](contributing.md#pull-request-workflow) to submit your 
changes for review.


### Editing the API Documentation

The API documentation is built from the docstrings in the Python `*.py` files under
the `pyatoa/pyatoa/` folder. **All docstrings** should follow the
[NumPy style guide](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard).
All functions/classes/methods should have docstrings with a full description of all
arguments and return values.

While the maximum line length for code is automatically set by Black, docstrings
must be formatted manually. To play nicely with Jupyter and IPython, **keep docstrings
limited to 79 characters** per line.


### Cross-referencing with Sphinx

The API reference is manually assembled in `doc/api/index.rst`.
The *autodoc* sphinx extension will automatically create pages for each
function/class/module/method listed there.

You can reference functions, classes, modules, and methods from anywhere
(including docstrings) using:

- <code>:func:\`Function name <package.module.function\>`</code>
- <code>:class:\`Class name <package.module.class>\`</code>
- <code>:meth:\`Method name <package.module.method>\`</code>
- <code>:mod:\`Module name <package.module>\`</code>

Sphinx will create a link to the automatically generated page for that
function/class/module/method.

## Contributing Code

### Pyatoa Code Overview

The source code for Pyatoa is located in the `pyatoa/` directory. When contributing
code, be sure to follow the general guidelines in the
[pull request workflow](contributing.md#pull-request-workflow) section.


### Testing your Code

Automated testing helps ensure that our code is as free of bugs as it can be.
It also lets us know immediately if a change we make breaks any other part of the code.

All of our test code and data are stored in the `tests` directory.
We use the [pytest](https://pytest.org/) framework to run the test suite.

Please write tests for your code so that we can be sure that it won't break any of the
existing functionality.
Tests also help us be confident that we won't break your code in the future.

If you're **new to testing**, see existing test files for examples of things to do.
**Don't let the tests keep you from submitting your contribution!**
If you're not sure how to do this or are having trouble, submit your pull request
anyway.
We will help you create the tests and sort out any kind of problem during code review.

You can also run tests in just one test script using:

    cd pyatoa/tests
    pytest NAME_OF_TEST_FILE.py

or run tests which contain names that match a specific keyword expression:

    cd pyatoa/tests
    pytest -k KEYWORD NAME_OF_TEST_FILE.py
