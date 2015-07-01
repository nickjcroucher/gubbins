require 'formula'

class Gubbins < Formula
  # tag "bioinformatics"
  # doi "10.1093/nar/gku1196"

  homepage 'https://github.com/andrewjpage/gubbins'
  url 'https://github.com/andrewjpage/gubbins/archive/convert_python_to_3.tar.gz'
  sha1 'e779b4f9def6e671587248b92c00fc8e60107314'
  head 'https://github.com/andrewjpage/gubbins.git'

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
