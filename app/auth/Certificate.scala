package auth

import java.io.FileInputStream
import java.security.cert.{Certificate, CertificateFactory}


private object Certificate {

  def read(filename: String): Certificate = {
    val fin = new FileInputStream(filename)
    val f = CertificateFactory.getInstance("X.509")
    f.generateCertificate(fin)
  }

}