package auth

import java.io.{FileInputStream, FileOutputStream, InputStreamReader, OutputStreamWriter}
import java.security._
import java.security.interfaces.{RSAPrivateKey, RSAPublicKey}
import java.security.spec.{PKCS8EncodedKeySpec, X509EncodedKeySpec}

import org.bouncycastle.jce.provider.BouncyCastleProvider
import org.bouncycastle.util.io.pem.{PemObject, PemReader, PemWriter}


// Ref:
// https://www.txedo.com/blog/java-generate-rsa-keys-write-pem-file/
// https://www.txedo.com/blog/java-read-rsa-keys-pem-file/

object RSA {

  Security.addProvider(new BouncyCastleProvider())
  val generator: KeyPairGenerator = KeyPairGenerator.getInstance("RSA", "BC")

  private def generateRSAKeyPair(keySize: Int = 1024): KeyPair = {
    generator.initialize(keySize)
    generator.generateKeyPair()
  }

  class RSAPair(val keySize: Int = 1024) {
    val keyPair: KeyPair = generateRSAKeyPair(keySize)
    val privateKey: RSAPrivateKey = keyPair.getPrivate.asInstanceOf[RSAPrivateKey]
    val publicKey: PublicKey = keyPair.getPublic.asInstanceOf[RSAPublicKey]
  }

  def writeTestKeyPair(keySize: Int = 1024, path: String = "resources/rsa_keys/test"): Unit = {
    val keyPair = new RSAPair(keySize)
    writePublicKeyToPemFile(keyPair.publicKey, path +"/id_rsa.pub")
    writePrivateKeyToPemFile(keyPair.privateKey, path +"/id_rsa")
  }

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

  def writePublicKeyToPemFile(key: PublicKey, path: String = "id_rsa.pub"): Unit = {
    PEMFile.write(key, "RSA PUBLIC KEY", path)
  }

  def writePrivateKeyToPemFile(key: PrivateKey, path: String = "id_rsa"): Unit = {
    PEMFile.write(key, "RSA PRIVATE KEY", path)
  }

  object PEMFile {

    def write(key: Key, description: String, filename: String): Unit = {
      val pemObject = new PemObject(description, key.getEncoded)
      val pemWriter: PemWriter = new PemWriter(new OutputStreamWriter(new FileOutputStream(filename)))
      try {
        pemWriter.writeObject(pemObject)
      } finally {
        pemWriter.close()
      }
    }

    def read(path: String): PemObject = {
      val pemReader: PemReader = new PemReader(new InputStreamReader(new FileInputStream(path)))
      try {
        val pemObject: PemObject = pemReader.readPemObject()
        pemObject
      } finally {
        pemReader.close()
      }
    }

  }

}
