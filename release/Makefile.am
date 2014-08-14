# Aidan Delaney <a.j.delaney@brighton.ac.uk> is responsible for this,
# blame him for any errors.

# Requires the use of Vagrant and reasonable bandwidth/time.

PACKAGE_VERSION=$(VERSION)~trusty2

release-local: ubuntu-release

# copy gubbins using /vagrant default mount rather than puppet, to avail of automake substitution.
ubuntu-release:
	vagrant up
	vagrant ssh -c "sudo apt-get update"
	vagrant reload
	vagrant provision
	cp ../${PACKAGE}-${VERSION}.tar.gz .
	vagrant ssh -c "tar xzvf /vagrant/${PACKAGE}-${VERSION}.tar.gz"
	vagrant ssh -c "cd ${PACKAGE}-${VERSION} && dpkg-buildpackage -uc -us -rfakeroot"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}_amd64.deb /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}_amd64.changes /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}.tar.gz /vagrant"
	vagrant halt

source:
	vagrant up
	vagrant provision
	cp ../${PACKAGE}-${VERSION}.tar.gz .
	vagrant ssh -c "tar xzvf /vagrant/${PACKAGE}-${VERSION}.tar.gz"
	vagrant ssh -c "mv ${PACKAGE}-${VERSION} ${PACKAGE}-${PACKAGE_VERSION}"
	vagrant ssh -c "cd ${PACKAGE}-${PACKAGE_VERSION} && dpkg-buildpackage -uc -us -rfakeroot"
	vagrant ssh -c "cd ${PACKAGE}-${PACKAGE_VERSION} && make clean"
	vagrant ssh -c "dpkg-source -rfakeroot -b ${PACKAGE}-${PACKAGE_VERSION}"
	vagrant ssh -c "cd ${PACKAGE}-${PACKAGE_VERSION} && dpkg-genchanges -S > ../${PACKAGE}_${PACKAGE_VERSION}_amd64.changes"
	vagrant ssh -c "cp ${PACKAGE}_${PACKAGE_VERSION}.dsc /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${PACKAGE_VERSION}_amd64.changes /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${PACKAGE_VERSION}.tar.gz /vagrant"
	vagrant halt