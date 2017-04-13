package forms

import models.Sample
import play.api.libs.json.JsValue


object SamplesForm {
  def fromJson(j: JsValue): Sample = {
    Sample(
      (j \ "name").as[String],
      (j \ "filename").as[String],
      (j \ "project").asOpt[String],
      (j \ "hash").asOpt[String],
      (j \ "description").asOpt[String],
      (j \ "isActive").asOpt[Boolean],
      (j \ "isAdmin").asOpt[Boolean] getOrElse false
    )
  }
}
