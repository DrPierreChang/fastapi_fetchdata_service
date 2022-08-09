from datetime import datetime
from typing import Union
import uuid as uuid_pkg

from pydantic import BaseModel, EmailStr, constr

__all__ = (
    "UserModel",
    "UserCreate",
    "UserLogin",
    "UserUpdate",
    "UserDetail",
)


class UserBase(BaseModel):
    username: constr(min_length=5, max_length=20, regex="^[a-zA-Z0-9_-]+$")


class UserLogin(UserBase):
    password: constr(min_length=5, max_length=20,)


class UserCreate(UserLogin):
    email: Union[EmailStr, None] = None


class UserModel(UserBase):
    uuid: uuid_pkg.UUID
    created_at: datetime
    is_superuser: bool
    is_active: bool
    is_totp_enabled: bool
    hashed_password: str
    email: Union[EmailStr, None] = None


class UserUpdate(BaseModel):
    username: Union[constr(min_length=5, max_length=20,
                           regex="^[a-zA-Z0-9_-]+$"), None] = None
    password: Union[constr(min_length=5, max_length=20,), None] = None
    email: Union[EmailStr, None] = None


class UserDetail(UserBase):
    uuid: uuid_pkg.UUID
    created_at: datetime
    is_superuser: bool
    is_active: bool
    is_totp_enabled: bool
    email: Union[EmailStr, None] = None
