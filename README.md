# Getting Started

It's recommended that you start by downloading the Anaconda Python Distribution (we will use Python 3 syntax for this class) from https://docs.conda.io/projects/conda/en/latest/user-guide/install/; it's available for Windows, Mac, and Linux, and it should come with all the packages you'll need for this class pre-installed.

Windows users may want to consider setting up VirtualBox with Ubuntu Linux. While this is not necessary for this class, scientific computing is significantly complicated by using a Windows machine, and you'll need to use a workaround to get the _git_ tool used for version control.

If you don't already have a Python IDE that you like, you can download the PyCharm IDE from https://www.jetbrains.com/pycharm/download/; it has a lot of features and is much  simpler to use than a text editor. You should be able to get free access to the Professional edition with your .edu email; if you don't have one, the Community edition is also sufficient.

If you have PyCharm installed, you'll want to go into Preferences and make sure that the anaconda python is set as the project interpreter for this project. Type "interpreter" into the searchbox and the project interpreter option should show up.

Finally, you'll want to set up git on your machine. If you have a Mac, type "git" at the command line, and you'll download some XCode Developer Tools that will set up git automatically; 

If you have a Windows machine, things might be a little trickier. I recommend just downloading git for windows from https://git-scm.com/download/win. This will come with it's own command line tool that will allow you to clone the repository using `git clone https://github.com/eeskin/CM122_starter_code.git`. 

# Installing dependencies

We recommend creating a new python environment for this class. The starter code doesn't have that many dependencies and things should come pre installed with your Anaconda distribution, but it's general good practice to keep things in separate environments, especially when you're running somebody else's code. Unix users can do this via command line by running
```
conda create -n CM122 python=3.7
conda activate CM122
```

You can then use `conda install package_name` to install any package you might need. The description of each project will mention if you need to install any packages for that project.

For Windows users and those who are generally adverse to using the command-line, you don't __really__ need to have a separate environment for this project, but if you wish to do it, you can manage virtual environments using PyCharm as well.
