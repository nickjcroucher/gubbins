require 'formula'

class Fastml2 < Formula
  homepage 'https://github.com/sanger-pathogens/fastml'
  url 'https://github.com/sanger-pathogens/fastml/archive/v2.3.tar.gz'
  sha1 'e70b2c83a321ae1664e1ab2a1d16b6bff70fdb01'
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

