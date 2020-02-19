# Contributing to GUDHI

First of all, thank you for the time you may take to contribute to GUDHI !

# In case you have a question

Please, check our [contact web page](https://gudhi.inria.fr/contact/).

# In case you found an issue

Please, first check [opened issues on GUDHI](https://github.com/GUDHI/gudhi-devel/issues).

If the problem you are facing is not referenced, do not hesitate to open a [new issue](https://github.com/GUDHI/gudhi-devel/issues/new).

This place is also a good place if you have some enhancement you want to propose.
There is a label **enhancement** in the [new issue](https://github.com/GUDHI/gudhi-devel/issues/new) page.

# In case you want to contribute to GUDHI

## You are not familiar with GitHub ?

Please take some time to read our [how to use GitHub to contribute to GUDHI](/home/vincent/workspace/gudhi/gudhi-devel/for_dev/how_to_use_github_to_contribute_to_gudhi.md).

## Something you want to improve in the documentation

For C++ documentation, you can find it in the directories:
* *src/common/doc* for the main page and installation instructions
* *src/NAME_OF_THE_MODULE/doc* for the main page of a module
* *src/NAME_OF_THE_MODULE/include/gudhi* for the documentation generated from the code.
We use Doxygen to generate the code and you will be able to verify the result in CircleCI Doxygen target in the artifacts.

For Python documentation, you can find it in the directories:
* *src/python/doc* for the main page, installation instructionsand for the main pages of the modules
* *src/python/gudhi/NAME_OF_THE_MODULE.pyx* for the documentation generated from the code.
We use Sphinx to generate the code and you will be able to verify the result in CircleCI Sphinx target in the artifacts.

## Something you want to improve in the code

We don't ask for any paperwork but we expect you don't submit anything you are not allowed to:
* check that your work contract and your employer allow you to contribute to this open source project.
* insure you do not violate someone's intellectual property.
* ...

Please, take some time to read our [code conventions](code_conventions.md)

As a convention, we set a Pull Request as a **Draft Pull Request** when we work on something we want the other contributors to see.

We click on **Ready for review** to ask for a peer review of the contribution.
