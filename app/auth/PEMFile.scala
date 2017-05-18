package auth

import java.io._
import java.security.Key
import org.bouncycastle.util.io.pem.{PemObject, PemReader, PemWriter}
import play.api.db.Database


object PEMFile {

  def read(path: String): PemObject = {
    val pemReader: PemReader = new PemReader(new InputStreamReader(new FileInputStream(path)))
    try {
      pemReader.readPemObject()
    } finally {
      pemReader.close()
    }
  }

  def write(key: Key, description: String, filename: String): Unit = {
    val pemObject = new PemObject(description, key.getEncoded)
    val pemWriter: PemWriter = new PemWriter(new OutputStreamWriter(new FileOutputStream(filename)))
    try {
      pemWriter.writeObject(pemObject)
    } finally {
      pemWriter.close()
    }
  }

  def writeToDb(appName: String, db: Database, key: Key, description: String): Unit = {
    val pemObject = new PemObject(description, key.getEncoded)
    db.withConnection { conn =>
      val stringWriter = new StringWriter()
      val pemWriter: PemWriter = new PemWriter(stringWriter)
      try {
        pemWriter.writeObject(pemObject)
      } finally {
        pemWriter.close()
      }
      val statement = conn.prepareStatement("""
          INSERT INTO `apps`(`keyFile`) VALUES (?) WHERE app = ?
        """)
      statement.setString(1, stringWriter.toString)
      statement.setString(2, appName)
    }
  }

}