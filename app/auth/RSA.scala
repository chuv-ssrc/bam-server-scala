package auth

import java.io._
import java.security._
import java.security.interfaces.{RSAPrivateKey, RSAPublicKey}
import java.security.spec.{PKCS8EncodedKeySpec, X509EncodedKeySpec}

import org.bouncycastle.jce.provider.BouncyCastleProvider
import org.bouncycastle.util.io.pem.PemObject
import play.api.db.Database

// Ref:
// https://www.txedo.com/blog/java-generate-rsa-keys-write-pem-file/
// https://www.txedo.com/blog/java-read-rsa-keys-pem-file/
// http://stackoverflow.com/questions/6358555/obtaining-public-key-from-certificate

object RSA {

  /******** PEM stuff *********/

  def readPublicKeyFromPemFile(filename: String): PublicKey = {
    val pemObject: PemObject = PEMFile.read(filename)
    val factory: KeyFactory = KeyFactory.getInstance("RSA", "BC")
    val content: Array[Byte] = pemObject.getContent
    val pubKeySpec: X509EncodedKeySpec = new X509EncodedKeySpec(content)
    factory.generatePublic(pubKeySpec)
  }

  def readPrivateKeyFromPemFile(path: String): PrivateKey = {
    val pemObject: PemObject = PEMFile.read(path)
    val factory: KeyFactory = KeyFactory.getInstance("RSA", "BC")
    val content: Array[Byte] = pemObject.getContent
    val privKeySpec: PKCS8EncodedKeySpec = new PKCS8EncodedKeySpec(content)
    factory.generatePrivate(privKeySpec)
  }

  def writePublicKeyToDb(appName: String, db: Database, key: PublicKey): Unit = {
    PEMFile.writeToDb(appName, db, key, "RSA PUBLIC KEY")
  }

  def writePublicKeyToPemFile(key: PublicKey, path: String = "id_rsa.pub"): Unit = {
    PEMFile.write(key, "RSA PUBLIC KEY", path)
  }

  def writePrivateKeyToPemFile(key: PrivateKey, path: String = "id_rsa"): Unit = {
    PEMFile.write(key, "RSA PRIVATE KEY", path)
  }


  /******** CRT stuff *********/

  def readPublicKeyFromCertificate(filename: String): PublicKey = {
    Certificate.read(filename).getPublicKey
  }

  
  /***** Generate test keys *****/

  Security.addProvider(new BouncyCastleProvider())
  val generator: KeyPairGenerator = KeyPairGenerator.getInstance("RSA", "BC")

  private def generateRSAKeyPair(keySize: Int): KeyPair = {
    generator.initialize(keySize)
    generator.generateKeyPair()
  }

  private class RSAPair(val keySize: Int = 1024) {
    val keyPair: KeyPair = generateRSAKeyPair(keySize)
    val privateKey: RSAPrivateKey = keyPair.getPrivate.asInstanceOf[RSAPrivateKey]
    val publicKey: PublicKey = keyPair.getPublic.asInstanceOf[RSAPublicKey]
  }

  def writetestSample1Pair(keySize: Int = 1024, path: String = "resources/rsa_keys/test"): Unit = {
    val keyPair = new RSAPair(keySize)
    writePublicKeyToPemFile(keyPair.publicKey, path +"/id_rsa.pub")
    writePrivateKeyToPemFile(keyPair.privateKey, path +"/id_rsa")
  }

}
