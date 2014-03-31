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