# Aidan Delaney <a.j.delaney@brighton.ac.uk> is responsible for this,
# blame him for any errors.

# Requires the use of Vagrant and reasonable bandwidth/time.

release-local: ubuntu-release

# copy gubbins using /vagrant default mount rather than puppet, to avail of automake substitution.
ubuntu-release:
	vagrant up
	vagrant reload
	vagrant provision
	cp ../${PACKAGE}-${VERSION}.tar.gz .
	vagrant ssh -c "tar xzvf /vagrant/${PACKAGE}-${VERSION}.tar.gz"
	vagrant ssh -c "cd ${PACKAGE}-${VERSION} && dpkg-buildpackage -uc -us -rfakeroot"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}_amd64.deb /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}_amd64.changes /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}.tar.gz /vagrant"
	vagrant halt
