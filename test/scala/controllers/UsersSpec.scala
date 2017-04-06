package scala.controllers

import org.scalatestplus.play._
import play.api.db.Database
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._
import scala.setup.{WithDatabase, WithToken}
import utils.DbUtils._


/**
  * Test IndexController.
  */
class UsersSpec extends PlaySpec with OneAppPerSuite with WithToken with WithDatabase {


  "Add a user" should {

    "add a row to the users table" in {
      val body = Json.parse("""
        { "users": [{"username": "A"}, {"username": "B"}] }
      """)
      val db = dbContext.db
      val count0 = countRows(db, "users")
      val response = route(app, FakeAuthorizedRequest(PUT, "/users").withJsonBody(body)).get
      status(response) mustBe OK
      val count = countRows(db, "users")
      count must equal(count0 + 2)
    }

  }

}
