package scala.setup

import play.api.Application
import play.api.db.{Database, Databases}
import javax.inject._


@Singleton
class DbContext @Inject()(val db: Database)


trait WithDatabase {

  protected def dbContext(implicit app: Application): DbContext = {
    Application.instanceCache[DbContext].apply(app)
  }

}
