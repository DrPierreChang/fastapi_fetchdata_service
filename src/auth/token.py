from datetime import datetime, timedelta
import jwt
import uuid as uuid_pkg

from fastapi import HTTPException, status

from src.core import config


class Token:
    @staticmethod
    def encode_token(uuid: str, token_type: str, expire: int):
        payload = {
            'exp': datetime.utcnow() + timedelta(seconds=expire),
            'iat': datetime.utcnow(),
            'type': token_type,
            'sub': uuid,
            'jti': str(uuid_pkg.uuid4())
        }
        return jwt.encode(
            payload,
            config.JWT_SECRET_KEY,
            algorithm=config.JWT_ALGORITHM
        )

    @staticmethod
    def decode_token(token: str):
        try:
            payload = jwt.decode(token,
                                 config.JWT_SECRET_KEY,
                                 algorithms=['HS256'])
            return payload
        except jwt.ExpiredSignatureError:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                                detail='Token expired',
                                headers={"WWW-Authenticate": "Bearer"}
                                )
        except jwt.InvalidTokenError:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                                detail='Invalid token',
                                headers={"WWW-Authenticate": "Bearer"}
                                )

    @staticmethod
    def check_token_type(token: dict, token_type: str = 'access'):
        if token['type'] != token_type:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED,
                                detail='Wrong token type',
                                headers={"WWW-Authenticate": "Bearer"}
                                )
