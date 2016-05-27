package { "dh-make":
    ensure => "installed"
    }

package { ["gcc", "build-essential", "pkg-config","ntp","libtool"]:
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
package {"language-pack-en":
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

package {"python-biopython":
    ensure => "installed"
}

package {"python-reportlab":
  ensure => "installed"
}

package {"python-nose":
  ensure => "installed"
}

package {"python3-pip":
  ensure => "installed"
}


# The Debian/Ubuntu system biopython library has no egg-info associated with it
# https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=743927
# so setuptools will pull down an egg and needs Python.h available to make it
# all work.
package {"python-dev":
  ensure => "installed"
}

include apt

