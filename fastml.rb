require 'formula'

class Fastml < Formula
  homepage 'https://github.com/sanger-pathogens/fastml'
  url 'https://github.com/sanger-pathogens/fastml/archive/v2.0.3.tar.gz'
  sha1 '1b6afc4debfd52d6de09bc70e425d6871eb8e90c'
  head 'https://github.com/sanger-pathogens/fastml.git'

  depends_on :autoconf
  depends_on :automake
  depends_on :libtool

  def install
    system "make 2>/dev/null"
    bin.install 'programs/fastml/fastml'
  end
  
  test do
    system 'fastml 2>&1'
  end
end