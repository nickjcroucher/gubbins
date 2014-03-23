# Aidan Delaney <a.j.delaney@brighton.ac.uk> is responsible for this,
# blame him for any errors.

# Requires the use of Vagrant and reasonable bandwidth/time.

release: ubuntu-release

# copy gubbins using scp rather than puppet, to avail of automake substitution.
ubuntu-release:
	vagrant up
	vagrant provision
	sshpass -p "vagrant" scp -P 2222 ../${PACKAGE}-${VERSION}.tar.gz vagrant@localhost:~
	vagrant ssh -c "tar xzvf ${PACKAGE}-${VERSION}.tar.gz"
	vagrant ssh -c "cd ${PACKAGE}-${VERSION} && dpkg-buildpackage -rfakeroot"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}_amd64.deb /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}_amd64.changes /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}.dsc /vagrant"
	vagrant ssh -c "cp ${PACKAGE}_${VERSION}.tar.gz /vagrant"
	vagrant halt


.PHONY: release
