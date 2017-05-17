package utils

import java.io.File
import java.nio.file.Paths
import java.security.MessageDigest
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.{Await, Future}
import scala.concurrent.duration._
import scala.util.matching.Regex
import sys.process._



object Common {

  /**
    * Hash a string
    * @param s: string to hash.
    * @param method: hashing algorithm.
    */
  def hash(s: String, method: String = "SHA") = {
    MessageDigest.getInstance(method).digest((s).getBytes)
      .map(0xFF & _).map { "%02x".format(_) }.foldLeft(""){_ + _}
  }

  /**
    * Use the script onDisk.sh to check if the file exists,
    * because it is able to deal with long-term storage and mounted volumes like we have with /archive.
    */
  def isOnDisk(file: File, archive: Boolean = false): Boolean = {
    if (archive) {
      val isFound: Int = s"scripts/onDisk.sh ${file.toString}".!
      isFound == 0
    } else {
      file.exists
    }
  }

  /**
    * List all files in a directory, including subdirectories, recursively.
    */
  def listFilesTree(dir: File): Array[File] = {
    val these = dir.listFiles
    these ++ these.filter(_.isDirectory).flatMap(listFilesTree)
  }

  /**
    * Find all files with extension ".<extension>" in the given path.
    */
  def findInTree(path: String, extension: String): Array[File] = {
    val dir = Paths.get(path).toFile
    val regex = """.*\.""" + Regex.quote(extension) + """$"""
    listFilesTree(dir).filter(f => regex.r.findFirstIn(f.getName).isDefined)
  }

  def await[E](f: Future[E]): E = {
    Await.result(f, 1.second)
  }

}



