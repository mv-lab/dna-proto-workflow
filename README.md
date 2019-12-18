# DNA Snakemake Proto-Workflow

This repository is based on the [Snakemake tutorial](http://snakemake.readthedocs.io/en/latest/tutorial/welcome.html) and [PaneucalyptShortReads](https://github.com/kdmurray91/PaneucalyptShortReads)


 <img src="https://divingintogeneticsandgenomics.rbind.io/img/snakemake.png" alt="Girl in a jacket" width="200" height="150">  

----

## [Setup](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html)

The easiest way to setup these prerequisites is to use the **Miniconda** Python 3 distribution.

```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```

**Creating an environment with the required software**

> NOTE: conda's enviroment name in these examples is `snkmk`.

```
$ conda env create --name snkmk --file envs/condaenv.yml
```

**Activate the enviroment using:**

```
$ conda activate snkmk
```

Check this [conda cheatsheet](https://gist.github.com/mv-lab/62318ff0023bd626f1e05ed9c0155fd5) for more useful commands and tricks

<br>

----


## Reproduce results

Clone this repository into your machine:

```
git clone https://github.com/mv-lab/dna-proto-workflow.git
```
Check the new directory
```
cd dna-proto-workflow
ls
```

> NOTE: the ```output``` directory and ```metadata``` files are generated by the workflow.

get **dag.svg** and check rules (see ```--dag -npr```)

```
snakemake --dag -npr -j -1 | dot -Tsvg > dag.svg
eog dag.svg
```

**Run complete Workflow**

```
snakemake  -npr -j -1
```

### Check config and metadata

```utils/check_config.py``` allows you to read and print the config files:

```
from utils.check_config import *
config = readconfig('config.yml')
pconfig(config)
config['qc']
```


```utils/check_metadata.py``` allows you to print the metadata generated at the ```common``` rule:

```
python utils/check_metadata.py
```

### Data

Main directory content:

```
.
├── config.yml
├── envs
├── genomes_and_annotations
├── metadata
├── output
├── rules
├── samples
├── scripts
├── Snakefile
├── snpEff.config
├── utils
├── ...
```

The workflow uses the data from the following folders:

- **genomes_and_annotations**
- **samples**

You can download example data for testing the workflow. [click here to download](https://drive.google.com/drive/folders/1kpJsghU-jNTSKC9uEB9khos390lZNROr?usp=sharing)


<br>

-----

## How to contribute

<img align="right" width="300" src="https://github.com/firstcontributions/first-contributions/raw/master/assets/fork.png" alt="fork this repository" />

### Fork this repository

Fork this repository by clicking on the fork button on the top of this page.
This will create a copy of this repository in your GitHub account (not in your computer).


### Clone the repository

<br>

<img align="right" width="300" src="https://i.ibb.co/yVWsByF/Screenshot-from-2019-12-18-10-38-25.png" alt="clone this repository" />

Now **clone the forked repository** to your machine. 
Go to your GitHub account, open the forked repository, click on the clone button and then click the *copy to clipboard* icon. The url is going to be like: ```https://github.com/your-username/dna-proto-workflow.git``` where `your-username` is your GitHub username.

Open a terminal and run the following git command:

```
git clone https://github.com/your-username/dna-proto-workflow.git
```

<img align="right" width="300" src="https://github.com/firstcontributions/first-contributions/raw/master/assets/copy-to-clipboard.png" alt="copy URL to clipboard" />

Once you've cloned your fork, you can edit your local copy. However, if you want to contribute, you'll need to create a new branch as follows.

### Create a branch

Change to the repository directory on your computer (if you are not already there):
> NOTE: you can't change the name of the directory.

```
cd dna-proto-workflow
```

You can check your branches and what branch you're currently in, using the ```git branch``` command.

Now create a branch using the `git checkout` command:
```
git checkout -b new-branch-name
```

For example:
```
git checkout -b development
```

From this point, you are in the new branch
If something goes wrong, you can remove a branch using `git branch -d name-of-the-branch`. Or you can go back to the `master`branch using `git checkout master`.


### Make changes and commit

Once you modify something, if you go to the project directory and execute the command `git status`, you'll see there are changes.
Add those changes to your branch using the `git add` command:

```
git add .
or
git add name_of_the_file_you_modified
```

Now commit those changes using the `git commit` command:
```
git commit -m "write a message"
```

### Push changes to GitHub

Push your changes from your local copy (in your machine) to your Github (remote repository) using the command `git push`:
```
git push origin your-branch-name
```
replacing `your-branch-name` with the name of the branch you created earlier (eg `development`).


### Submit your changes for review

If you go to your repository on GitHub, you'll see a  `Compare & pull request` button. Click on that button.

<img style="float: right;" src="https://i.ibb.co/N7np2Ch/compare-and-pull.png" />

Now submit the pull request and you'll see something like:

<img style="float: right;" src="https://help.github.com/assets/images/help/pull_requests/pull-request-review-edit-branch.png" alt="submit pull request" />

We'll check your changes and merge them into this project (in general, into the `master` branch). You will get a notification email once the changes have been merged :)
