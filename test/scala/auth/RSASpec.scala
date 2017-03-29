package scala.auth

import org.scalatestplus.play._
import java.security._
import auth.RSA._



class RSASpec extends PlaySpec {

  "Read a private key from PEM file" in {
    val privateKey: PrivateKey = readPrivateKeyFromPemFile("resources/rsa_keys/test/id_rsa")
  }

  "Read a public key from PEM file" in {
    val publicKey: PublicKey = readPublicKeyFromPemFile("resources/rsa_keys/test/id_rsa.pub")
  }

  "Read a public key from Auth0 certificate" in {
    val auth0PublicKey: PublicKey = readPublicKeyFromCertificate("resources/rsa_keys/auth0.cer")
  }

}
