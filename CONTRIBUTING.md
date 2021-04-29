OSCILOS is being developed for research. Please [join us](http://www.oscilos.com/)!. 

## Git and Pull requests
* Contributions are welcome. For outside users, they can be submitted with GitHub pull request, and will be reviewed and accepted by the team. The project uses the Fork & Pull model; for more details see [here](https://help.github.com/articles/using-pull-requests).
* Internal developers have write access to the branches, and can merge changes themselves. They do not need to make Pull requests.
* The project uses the [GitFlow branching model](http://nvie.com/posts/a-successful-git-branching-model/). As such please make a new branch for every feature you are working on.
* Internal developers should merge Bugfixes and Hotfixes to `master-private` and `develop`.
* After a release is made to `master-private`, or a Bugfix/Hotfix is implemented, changes in `master-private` should be merged in the public repo. `master` branch.
* The latest changes available to users outside of MorgansLab are in the public repo. `master` branch. To have access to newer bug fixes or become a collaborator, contact the repo. admin. 
* Make sure your commit messages are clear. The first line should be a summary of the commit, following lines should contain details about the mentioned changes. 
* Please avoid using re-basing and fast forward merges. 

## Coding style and version numbering
* Use a lot of comments; more is better than not enough. 
* Indent your code using the Matlab automated indent feature. 
* When incrementing the version number, please stick to the [Semantic Versioning Standards](http://semver.org/)

## Documentation
* The documentation can be found in the [docs](docs) folder.
