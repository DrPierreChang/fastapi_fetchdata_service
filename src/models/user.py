from datetime import datetime
import uuid as uuid_pkg

from pydantic import EmailStr
from sqlalchemy import UniqueConstraint
from sqlmodel import Field, SQLModel

__all__ = ("User",)


class User(SQLModel, table=True):
    __table_args__ = (UniqueConstraint("username"),)
    uuid: uuid_pkg.UUID = Field(
        default_factory=uuid_pkg.uuid4,
        primary_key=True,
        index=True,
        nullable=False,
    )
    username: str = Field(nullable=False, min_length=5, max_length=20,
                          regex="^[a-zA-Z0-9_-]+$")
    email: EmailStr = Field(nullable=False)
    created_at: datetime = Field(default=datetime.utcnow(), nullable=False)
    is_superuser: bool = Field(default=False)
    is_active: bool = Field(default=True)
    is_totp_enabled: bool = Field(default=False)
    hashed_password: str = Field(nullable=False)
