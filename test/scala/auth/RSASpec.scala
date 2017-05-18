package scala.auth

import org.scalatestplus.play._
import java.security._

import auth.RSA._

import scala.setup.WithDatabase



class RSASpec extends PlaySpec with OneAppPerSuite with WithDatabase {

  val db = dbContext.db

  "Read a private key from PEM file" in {
    val privateKey: PrivateKey = readPrivateKeyFromPemFile("resources/rsa_keys/test/id_rsa")
  }

  "Read a public key from PEM file" in {
    val publicKey: PublicKey = readPublicKeyFromPemFile("resources/rsa_keys/test/id_rsa.pub")
  }

  "Read a public key from Auth0 certificate" in {
    val auth0PublicKey: PublicKey = readPublicKeyFromCertificate("resources/rsa_keys/auth0.cer")
  }

  //----- Storing the public key in db: ----//

  "Read a public key from a manual one-line string with returns" in {
    val keyString = "-----BEGIN RSA PUBLIC KEY-----\nMFwwDQYJKoZIhvcNAQEBBQADSwAwSAJBAK7ttYaE/1ldsb0OJQDQhhDWqwuFWIyt\nxgYIJH1HYA4UpA/Nm24fERIA1xi2Pomep6VTnQ/ThFP5hn2NyITwCIsCAwEAAQ==\n-----END RSA PUBLIC KEY-----"
    val publicKey: PublicKey = readPublicKeyFromString(keyString)
  }

  "Read a public key as string from PEM file" in {
    val source = scala.io.Source.fromFile("resources/rsa_keys/test/id_rsa.pub")
    val keyString = try source.getLines mkString "\n" finally source.close
    val publicKey: PublicKey = readPublicKeyFromString(keyString)
  }

  "Read a public key from database" in {
    val publicKey: PublicKey = readPublicKeyFromDb(db, 1)
  }

  "Write a public key as string to database" in {
    val keyString = "-----BEGIN RSA PUBLIC KEY-----\nMFwwDQYJKoZIhvcNAQEBBQADSwAwSAJBAK7ttYaE/1ldsb0OJQDQhhDWqwuFWIyt\nxgYIJH1HYA4UpA/Nm24fERIA1xi2Pomep6VTnQ/ThFP5hn2NyITwCIsCAwEAAQ==\n-----END RSA PUBLIC KEY-----"
    val publicKey: PublicKey = readPublicKeyFromString(keyString)
    writePublicKeyToDb(db, 1, publicKey)
  }

  "Write a public key from PEM file to database" in {
    val source = scala.io.Source.fromFile("resources/rsa_keys/test/id_rsa.pub")
    val keyString = try source.getLines mkString "\n" finally source.close
    val publicKey: PublicKey = readPublicKeyFromString(keyString)
    writePublicKeyToDb(db, 1, publicKey)
  }

  "Write a public key from certificate to database" in {
    val auth0PublicKey: PublicKey = readPublicKeyFromCertificate("resources/rsa_keys/auth0.cer")
    writePublicKeyToDb(db, 1, auth0PublicKey)
  }

}
