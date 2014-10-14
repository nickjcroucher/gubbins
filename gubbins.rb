require 'formula'

class Gubbins < Formula
  homepage 'https://github.com/sanger-pathogens/gubbins'
  url 'https://github.com/sanger-pathogens/gubbins/archive/v1.0.1.tar.gz'
  sha1 '0a4e75fb511de2270f4ad16390fb090d6ddb6ad7'
  head 'https://github.com/sanger-pathogens/gubbins.git'

  depends_on :autoconf
  depends_on :automake
  depends_on :libtool
  depends_on 'raxml'
  depends_on 'fasttree'
  depends_on 'fastml'
  depends_on :python
  
  depends_on LanguageModuleDependency.new :python, "biopython", "Bio"
  
  def install
    inreplace "src/Makefile.am", "-lrt", "" if OS.mac? # no librt for osx
    inreplace "configure.ac", "PKG_CHECK_MODULES([zlib], [zlib])", "AC_CHECK_LIB(zlib, zlib)" if OS.mac?
    inreplace "python/Makefile.am", "--root=$(DESTDIR)", "--root=#{prefix}" if OS.mac?
    
    system "cat python/Makefile.am"
    system "autoreconf -i"
    system "./configure",
           "--disable-debug",
           "--disable-dependency-tracking",
           "--prefix=#{prefix}"

    system "make","install"
  end

  test do
    system "gubbins"
    system "run_gubbins.py"
  end
end

