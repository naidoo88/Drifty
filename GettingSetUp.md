# Setting up a virtual environment and Jupyter
## Context for the Python Noob

### Python vs. Python3
To cut a very long story short, Python(2.7) has been depreciated, and all new code should be developed using Python3.  Due to many packages still having dependencies on old python libraries, you will commonly find both included with your operating system.

`which python && python --version`

`which python3 && python3 --version`

We will be using Python3.

### pip
pip is the Package Installer for Python.  We will use this to install all of the required python packages.  If it is not already installed on your system (`which pip3` -- the '3' specifying the python3-flavour of pip), you will need to install this. 

eg. `sudo apt install python3-pip`

If your operating system doesn't use `apt` search how to install this for your specific case.  

__*Remember*, you want pip for Python3!__


## Why a virtual environment?
In short, because it will save you a lot of grief down the line. 

In longer form, lots of different scripts and applications rely on Python.  Therefore, as a general rule, you want to leave your global installation alone.  Furthermore, Python libraries are actively being developed all the time.  Every now and again you may find your old scripts do not work with a new version of a package.  Virtual environments allow you to set up a 'sand-box' with specific versions of a selection of packages with just a few commands.

There are various solutions to virtual environments, but here I present what I believe is the simplest solution to get you rolling and making corrections.

## Setting up the virtual environment with pip
  1. Clone the repository if you haven't already, and `cd` into the directory.

  2. Create a virtual environment:   
     `python3 -m venv "my_env_name"`  
     Choose whatever name you like.  This will create a directory of that name where our virtual environment is stored.  
     Depending on your python3 installation you may need to install venv (eg. `sudo apt install python3-venv`).
  
  3. Activate the virtual environment:  
     `source my_env_name/bin/activate`  
     There will also be an `activate.csh` and `activate.fish` variant if you are not using bash.
  
  4. Install all of the required packages in a one-liner:  
     `pip install -r requirements.txt`  
     Don't panic if you get errors for some of the packages on the list.  The most likely explanation is that the version of Python3 you are using has now incorporated those packages into the base install, and everything will work as expected.


## Running Jupyter-Notebook
With the above steps completed, you should have everything you need to run the tutorial notebook! 

With your virtual environment active, run:  `jupyter-notebook`.  This should automatically take you to a browser tab, showing you the directory of where you launched jupyter-notebook from.  

Click on `simkincorr_tutorial.ipynb` and you're off to the races!

**NB:** The line: `%matplotlib notebook` in the first code cell (under the imports) is important! Without it, the plots will not interactively update when you tweak them.

## Can I use Jupyter remotely?
Assuming no barriers have been set up by your sys-admins: **Yes!**
  
   1. `ssh` onto the remote machine which you want to use and set up the virtual environment, as described above.

   2. When running the notebook, we will use a couple of extra flags:  
      `jupyter-notebook --no-browser --port=8887`  
      * Instead of launching a browser, you will be provided with a hyperlink that you can use to access the notebook.
      * `--port=` allows us to specify which port we want Jupyter to use.  
         If you are using a shared machine, then as a courtesy to others who may want to use   Jupyter, we can choose something other than the default (8888).  
   
   3. On your local machine: forward the port you used above:  
      `ssh -NL 8887:localhost:8887 user@remote-host`  
         * `N` indicates that we do not actually want to access the remote  
         * `L` is the command for port-forwarding

   4. On your local machine:  open your browser of choice and paste in the link you were provided.

## Optional "quality of life" tweak
The Jupyter-Notebook GUI has a preset maximum cell-width (to match A4 paper), which can leave loads of unused space in a maximised browser.  As we are dealing with single, large figures, this means you may need to scroll to view the full width of the figure.

This can be fixed really easily, in just a few steps:

   1. Close Jupyter if it is already running.

   2. In your terminal, execute:  
   `jupyter --config-dir`  
      This will return the path to where jupyter looks for custom settings.  (eg `./home/.jupyter`)
   
   3. Create a directory at this location called 'custom':  
     `mkdir ./home/.jupyter/custom`

   4. Create a file called `./home/.jupyter/custom/custom.css` which contains the following:  
      ```css
      /* Make the notebook cells take almost all available width */
      .container {
          width: 99% !important;
      }   

      /* Prevent the edit cell highlight box from getting clipped;
       * important so that it also works when cell is in edit mode*/
      div.cell.selected {
          border-left-width: 1px !important;
      }
      ```
   5. When you launch jupyter-notebook you will now have full width cells, which adjust to the size of your browser window! 

## "Help, I get errors when I run the first cell!"
I include this because I have come across troubles in my time when trying to run Jupyter on a machine with a "busy" command line environment (eg. on a PC on your institution's network).

First things first: double check you activated your virtual environment before you launched Jupyter! 

Second: Check the package *is* installed!
   1. Note which package is throwing the error (eg, uproot).
   2. In your terminal, with your virtual environment active: `pip show uproot`   
      If the package is properly installed, you should see something like this:   
      ```csh
      (my-env) [pauln@npc* CalcCorrections]$ pip show uproot
      Name: uproot
      Version: 3.13.0
      Summary: ROOT I/O in pure Python and Numpy.
      Home-page: https://github.com/scikit-hep/uproot
      Author: Jim Pivarski (IRIS-HEP)
      Author-email: pivarski@princeton.edu
      License: BSD 3-clause
      Location: ./my-env/lib/python3.6/site-packages
      Requires: numpy, uproot-methods, awkward, cachetools
      Required-by: 
      ```

If not, something appears to have not worked when you installed the packages with the `requirements.txt` - try revisiting this step above. 

If the package shows as installed, it is likely that your environment is messing with you. 


Simply running jupyter-notebook when your environment is active *should* cause Jupyter to take the environment as the default kernel.  The "kernel" defines what Jupyter uses to run your code when you run a cell.

If this is not happening, we can be more explicit.  We do this by adding the our virtual environment to Jupyter as a kernel which we can then select:
   1. Kill Jupyter if it is still running (the process, not the browser tab!) by hitting `Ctrl+C` twice in quick succession in the terminal from which you launched it.
   2. Make sure your virtual environment still active.
   3. Add the virtual environment to the list of installed Jupyter kernels:  
      `python -m ipykernel install --user --name my-env --display-name "Python3(my-env)"`  
      where `--name` is what you named your virtual environment, and `--display-name` is what it will be called in the Jupyter GUI.
   4. Restart Jupyter and open `simkincorr_tutorial.ipynb`
   5. From the menu at the top: Kernel >> Change kernel >> Python3(my-env)




