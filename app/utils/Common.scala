package utils

import java.io.File
import java.security.MessageDigest
import play.api.libs.json.JsValue
import play.api.mvc.{AnyContent, Request}
import sys.process._



object Common {

  /**
    * Get a String value from the body of a POST request, given its key.
    */
  def stringValueFromRequestBody(request: Request[AnyContent], key: String): Option[String] = {
    request.body.asJson flatMap { body: JsValue =>
       (body \ key).asOpt[String]
    }
  }

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

}



