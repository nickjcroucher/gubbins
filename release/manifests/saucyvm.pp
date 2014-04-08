package { "dh-make":
    ensure => "installed"
    }

package { "gcc":
    ensure => "installed"
    }

package {"autotools-dev":
    ensure => "installed"
    }

package {"zlib1g-dev":
    ensure => "installed"
    }

package {"check":
    ensure => "installed"
    }

package {"python-setuptools":
    ensure => "installed"
    }

package {"cdbs":
    ensure => "installed"
    }

package {"raxml":
    ensure => "installed"
    }

package {"fasttree":
    ensure => "installed"
    }

include apt

# we need to pull in a packaged version of fastml for building.
# Supplied by Aidan Delaney <aidan@ontologyengineering.org>, so blame him.
apt::ppa { 'ppa:a-j-delaney/gubbins-ppa': }

package {"fastml":
        ensure => "installed"
        }
