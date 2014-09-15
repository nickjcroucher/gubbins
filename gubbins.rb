require 'formula'

class Gubbins < Formula
  homepage 'https://github.com/sanger-pathogens/gubbins'
  url 'https://github.com/sanger-pathogens/gubbins/archive/v0.6.tar.gz'
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

    system "autoreconf -i"
    system "./configure",
           "--disable-debug",
           "--disable-dependency-tracking",
           "--prefix=#{prefix}"

    system "make","DESTDIR=#{prefix}", "install"    
  end

  test do
    system "gubbins"
    system "run_gubbins.py"
  end
end