# Wait for Neo4j to be available
MAX_RETRIES=30
RETRY_COUNT=0

echo "Waiting for Neo4j to become available..."
until cypher-shell -u $NEO4J_USER -p $NEO4J_PASSWORD "RETURN 1" > /dev/null 2>&1 || [ $RETRY_COUNT -eq $MAX_RETRIES ]; do
  sleep 2
  RETRY_COUNT=$((RETRY_COUNT + 1))
  echo "Still waiting ($RETRY_COUNT/$MAX_RETRIES)..."
done

if [ $RETRY_COUNT -eq $MAX_RETRIES ]; then
  echo "Timed out waiting for Neo4j."
  exit 1
fi

echo "Creating database '$BC_TABLE_NAME'"
cypher-shell -u $NEO4J_USER -p $NEO4J_PASSWORD "CREATE DATABASE $BC_TABLE_NAME IF NOT EXISTS;"
echo "Database created or already exists!"
