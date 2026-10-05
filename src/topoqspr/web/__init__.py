"""Web interface on the v1.0 pipeline (exploratory mode; docs/WEB.md).

``topoqspr.web.app:app`` is the ASGI application (``uvicorn topoqspr.web.app:app``);
:mod:`topoqspr.web.jobs` holds the persistent job store and the background worker.
Needs the ``web`` extra (fastapi, uvicorn).
"""
