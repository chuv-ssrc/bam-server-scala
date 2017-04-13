package utils

import java.sql.ResultSet

import play.api.db.Database

import scala.concurrent.{Await, Future}
import scala.concurrent.duration._


object DbUtils {

  def countRows(db: Database, tableName: String): Int = {
    db.withConnection { conn =>
      val statement = conn.prepareStatement(s"SELECT COUNT(*) as count FROM `$tableName`;")
      val res = statement.executeQuery()
      res.next()
      res.getInt("count")
    }
  }

  def await[E](f: Future[E]): E = {
    Await.result(f, 1.second)
  }

}
