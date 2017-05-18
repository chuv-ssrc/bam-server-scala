package auth

import java.io._
import java.security.Key
import org.bouncycastle.util.io.pem.{PemObject, PemReader, PemWriter}
import play.api.db.Database


object PEMFile {

  def readFromFile(path: String): PemObject = {
    val pemReader: PemReader = new PemReader(new InputStreamReader(new FileInputStream(path)))
    try {
      pemReader.readPemObject()
    } finally {
      pemReader.close()
    }
  }

  def writeToFile(key: Key, description: String, filename: String): Unit = {
    val pemObject = new PemObject(description, key.getEncoded)
    val pemWriter: PemWriter = new PemWriter(new OutputStreamWriter(new FileOutputStream(filename)))
    try {
      pemWriter.writeObject(pemObject)
    } finally {
      pemWriter.close()
    }
  }

  def readFromString(keyString: String): PemObject = {
    val stringReader = new StringReader(keyString)
    val pemReader: PemReader = new PemReader(stringReader)
    try {
      pemReader.readPemObject()
    } finally {
      pemReader.close()
    }
  }

  def writeToString(key: Key, description: String): String = {
    val pemObject = new PemObject(description, key.getEncoded)
    val stringWriter = new StringWriter()
    val pemWriter: PemWriter = new PemWriter(stringWriter)
    try {
      pemWriter.writeObject(pemObject)
    } finally {
      pemWriter.close()
    }
    stringWriter.toString
  }

  def readFromDb(db: Database, appId: Int): PemObject = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement("""
          SELECT `key` FROM `apps` WHERE id = ?;
        """)
      statement.setInt(1, appId)
      val res = statement.executeQuery()
      res.next()
      val keyString = res.getString("key").replace("\\n", "\n")
      readFromString(keyString)
    }
  }

  def writeToDb(appId: Int, db: Database, key: Key, description: String): Unit = {
    db.withConnection { conn =>
      val keyString = writeToString(key, description)
      val statement = conn.prepareStatement("""
          UPDATE `apps` SET `key` = ? WHERE id = ?;
        """)
      statement.setString(1, keyString)
      statement.setInt(2, appId)
      statement.execute()
    }
  }

}