require 'formula'

class Gubbins < Formula
  # tag "bioinformatics"
  # doi "10.1093/nar/gku1196"

  homepage 'https://github.com/sanger-pathogens/gubbins'
  url 'https://github.com/sanger-pathogens/gubbins/archive/v1.4.0.tar.gz'
  sha1 'adf2d25eb1740155370f3c33c07c1622d4703c0c'
  head 'https://github.com/sanger-pathogens/gubbins.git'

  depends_on :autoconf
  depends_on :automake
  depends_on :libtool
  depends_on 'check'
  depends_on 'raxml'
  depends_on 'fasttree'
  depends_on 'fastml' => :recommended
  depends_on :python3
  depends_on 'zlib'
  
  def install
    inreplace "src/Makefile.am", "-lrt", "" if OS.mac? # no librt for osx
    inreplace "configure.ac", "PKG_CHECK_MODULES([zlib], [zlib])", "AC_CHECK_LIB(zlib, zlib)" if OS.mac?
    
    system "cat python/Makefile.am"
    system "autoreconf -i"
    system "./configure",
           "--disable-debug",
           "--disable-dependency-tracking",
           "--prefix=#{prefix}"
    system "make","check"
    system "make","install"
  end

  test do
    system "gubbins"
    system "run_gubbins.py"
    system "gubbins_drawer.py"
  end
end
