package scala.auth

import org.scalatestplus.play._
import play.api.libs.json._
import scala.setup.WithDatabase
import pdi.jwt.{JwtAlgorithm, JwtJson}



class HMACSpec extends PlaySpec with OneAppPerSuite with WithDatabase {

  val db = dbContext.db

  "Verify a HS256 jwt as returned by pyjwt (e.g. Varapp)" in {
    // >>> token = jwt.encode({"info": "test"}, "secret", algorithm='HS256').decode('utf-8')
    // >>> payload = jwt.decode(token, secret, algorithms=['HS256'], verify=True)
    val pytoken = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpbmZvIjoidGVzdCJ9.NerbtlDCuVHemOMgOYOzSbqMmVe8hYJMWKshgV4v2QI"
    JwtJson.validate(pytoken, "secret", Seq(JwtAlgorithm.HS256))
  }

  "Verify a HS256 jwt as returned by JwtJson.encode" in {
    val token2 = JwtJson.encode(Json.obj("info" -> "test"), "secret", JwtAlgorithm.HS256)
    JwtJson.validate(token2, "secret", Seq(JwtAlgorithm.HS256))
  }

}
