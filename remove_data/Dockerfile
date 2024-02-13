LABEL sushi.pipline.version = '0.0.1'
LABEL sushi.pipeline.vendor = 'PRISM'

FROM python:3.11-slim

WORKDIR /
COPY . .
ENTRYPOINT [ "aws_batch" ]

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt