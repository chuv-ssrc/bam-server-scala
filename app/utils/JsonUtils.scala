package utils

import play.api.libs.json._

object JsonUtils {

  implicit val anyWrites = new Writes[Any] {
    def writes(x: Any) = x match {
      case n: Int => JsNumber(n)
      case s: String => JsString(s)
      case b: Boolean => JsBoolean(b)
      case m: Map[String, Any] @unchecked => JsObject(m.mapValues(writes))
      case a: Seq[Any] => JsArray(a.map(writes))
      case _ => JsNull
    }
  }

}
