[![impacts.wiki](./PetaviusLangrenus_Poupeau_3000.png)](https://impacts.wiki)
# pyKO hydrocode

pyKO is a one-dimensional Lagrangian elastic-plastic hydrocode written in python. The code uses the Von Neumann finite difference scheme with second order accuracy. 
This python version of the KO code was developed from the book <a href="https://link.springer.com/book/10.1007/978-3-662-03885-7">Computer Simulation of Dynamic Phenomena by Mark Wilkins</a> (Springer-Verlag, 1999) and the <a href="https://www.eng.mu.edu/shockphysics/KO/">fortran KO code v11 by John Borg</a>.

This resource is available as part of the <a href="https://impacts.wiki">Impacts Community Wiki Project</a> under the GNU General Public License v3.0. This code is an easily modifiable and expandable teaching/learning tool that can be run in a Jupyter notebook or python interpreter.

The code is currently in a beta public release for testing and feedback at https://github.com/ImpactsWiki/pyko. 

## Help and Mailing List
If you use this code (in any programming language), join the user mailing list and post questions there:
ko-code-users@ucdavis.edu

To subscribe
* Send email to sympa@ucdavis.edu from the email address you want to subscribe
* Subject line: subscribe ko-code-users
* Body of email: leave empty
* You will receive an email providing a link to confirm your subscription

To unsubscribe
* Send email to sympa@ucdavis.edu from the email address to unsubscribe
* Subject: unsubscribe ko-code-users
* Body of email: leave empty

Report bugs in the python code on GitHub at https://github.com/ImpactsWiki/pyko/issues

## Other versions of the KO hydrocode

The KO code is available in multiple programming languages, with different features implemented in each version:
* <a href="https://www.eng.mu.edu/shockphysics/Workshops/">Fortran version (currently at v13) by John Borg</a>. v11 is included in the pyKO repository as it is used in the testing/tutorial examples.
* <a href="https://github.com/Yatagarasu50469/KO-Hydrocode">C version by David Helminiak</a>
* Matlab version by <a href="https://www.westpoint.edu/civil-and-mechanical-engineering/profile/nathaniel_helminiak-phd-eit">Nathan Helminiak</a> is included in the pyKO repository. 
