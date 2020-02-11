# How to use github to contribute to gudhi

Similar information is available in many places:
* https://jarv.is/notes/how-to-pull-request-fork-github/ (this one is using `upstream/master` when creating a new branch)
* https://help.github.com/en/github/getting-started-with-github/fork-a-repo
* https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/
* https://gist.github.com/Chaser324/ce0505fbed06b947d962
* etc

## Get a github account
I assume the account is called **LOGIN**, please replace as appropriate below. Log in to github.com using this account.

## Fork GUDHI/gudhi-devel project
Go to https://github.com/GUDHI/gudhi-devel and click on **fork** (top right).
Feel free to also click on the star next to it to show you like the project!
You can see your fork at https://github.com/LOGIN/gudhi-devel

## Create a local clone on your computer
```bash
git clone https://github.com/LOGIN/gudhi-devel.git
```

This creates a directory gudhi-devel, which you are free to move around or rename. For the following, change to that directory:
```bash
cd gudhi-devel
```

## Configuring a remote for a fork
```bash
git remote add upstream https://github.com/GUDHI/gudhi-devel.git
```

because you want to see the real gudhi, not just your clone.
(It is perfectly possible to do things in the reverse order, clone from GUDHI and add the one in LOGIN as extra remote, but the names of the remotes may not match the rest of this document. You can change the name of a remote with `git remote rename oldname newname`)

## Optional remotes
Optional, if you are interested in one of the old branches
```bash
git remote add oldies https://github.com/GUDHI/branches.git
```

Or if you want to spy on someone's work. I assume the someone's account is called **SOMEONE**
```bash
git remote add someone https://github.com/SOMEONE/gudhi-devel.git
```

## Stay up-to-date
```bash
git fetch -p --all
```
This is a command you can run quite regularly.
It tells git to check all that happened on github.
It is safe, it will not mess with your files.

## Create a branch, based on the current master
```bash
git checkout -b some-fancy-name --no-track upstream/master
```
Your local branch `master` and the one on your github clone are useless and often outdated, but for technical reasons there has to exist at least one branch at all times, it might as well be that one. upstream/master is the real deal, that's what you want to base your new branch on.

## The real coding is here!
Edit files, test, etc.

## Commit your changes (locally)
The basic command is just `git commit`, but it will do nothing by default.
You need `git add my_new_file` for every new file you want to commit.
And usually you'll want to use `git commit -a` so that all files that git already knows about and that have been modified get committed.

## Push your changes (remotely)
```bash
git push -u origin some-fancy-name
```
This puts a copy of your branch on your online clone of gudhi-devel.
Because of `-u`, it will remember where you like to push this branch, and next time you can just use `git push`.

## Play again!
Possibly iterate a few times, add more commits and push them.

## Your pull request is ready
Get your web browser to https://github.com/LOGIN/gudhi-devel, click on the button that says **Branch: some-name** (below the number of commits, above the list of files) and select the branch you are so proud of.
Click on **New pull request** next to it.

## Follow the instructions ;-)
Note that if your branch is not quite ready, you can make a **draft pull request** (see the arrow next to the confirmation button), and later you will have access to a button to say that the branch is ready for reviews now.
Draft pull requests can be a way to advertise that you are working on something, and possibly ask others for comments or help.

## Code review
Make sure you follow the discussion on your pull request, answer questions, take comments into account.
You can keep pushing new commits on your branch to your fork of gudhi-devel, the pull request will automatically notice the new commits there.
There is no need to create a new pull request.
Once the branch is under review, fixing issues is good, but please refrain from adding extra features, that just makes the reviewers' job harder and thus slower.
You may want to look at https://github.com/settings/notifications (and other settings nearby) if you don't receive emails when people comment on your pull request.
Some bold reviewer might make changes to your branch. You will then need `git pull` for your local branch to reflect those.

## Your work is merged!
Once your pull request has been closed (your branch merged), you can remove your branch, both locally and also the branch on your github fork:
```bash
git checkout master # or any other branch, but you cannot remove the branch you are currently in
git branch -d some-fancy-name # local branch delete
git push origin --delete some-fancy-name # remote branch delete
```
If you add @VincentRouvreau or @mglisse as collaborator (https://github.com/LOGIN/gudhi-devel/settings/collaboration), they may remove the branch on your clone at the same time as they merge the branch, so you only have the local one to remove (or keep if you are nostalgic).

## Keep in touch
Create a new branch and keep contributing!

Do not try to reuse an old branch that has already been merged.
Make sure you run the fetch command just before creating any new branch, so you don't base it on some outdated version of master.
You can also work on several branches at the same time, using `git checkout some-fancy-name` and `git checkout name-of-other-branch` to switch between them (commit before switching or things may get complicated).
